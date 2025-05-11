#!/bin/bash

while true; do
  echo "========================================"
  echo " Git 操作選單"
  echo "========================================"
  echo "1) 初始化主專案 remote 與 submodules（已 clone）"
  echo "2) 初始化 submodules（含檢查）"
  echo "3) 從 upstream 更新（含 submodules）"
  echo "4) 推送至 GitHub（含 submodules）"
  echo "X) 離開"
  echo "========================================"
  read -p "請輸入操作項目: " choice

  case "$choice" in
    1)
      echo "[INFO] 設定主 repo upstream..."
      if ! git remote get-url upstream > /dev/null 2>&1; then
        read -p "請輸入 upstream URL: " upstreamURL
        git remote add upstream "$upstreamURL"
      else
        echo "[INFO] upstream 已存在，略過。"
      fi

      if [ -f "submodules_config.txt" ]; then
        echo "[INFO] 初始化 submodules..."
        while read -r line; do
          set -- $line
          subPath="$1"
          originURL="$2"
          upstreamURL="$3"
          if [ -z "$(git config -f .gitmodules --get submodule.${subPath}.url)" ]; then
            git submodule add "$originURL" "$subPath"
          else
            echo "[INFO] submodule $subPath 已存在，略過。"
          fi
        done < submodules_config.txt
        git submodule init
        git submodule update --remote --merge
      else
        echo "[INFO] 無 submodules_config.txt，略過。"
      fi
      ;;
    2)
      if [ -f ".gitmodules" ]; then
        read -p "[INFO] 偵測到 submodule 記錄，是否重新初始化？(Y/N): " redo
        if [[ ! "$redo" =~ ^[Yy]$ ]]; then
          echo "[INFO] 已取消初始化。"
          continue
        fi
      fi

      if [ -f "submodules_config.txt" ]; then
        echo "[INFO] 開始初始化 submodules..."
        while read -r line; do
          set -- $line
          subPath="$1"
          originURL="$2"
          upstreamURL="$3"
          if [ -z "$(git config -f .gitmodules --get submodule.${subPath}.url)" ]; then
            git submodule add "$originURL" "$subPath"
          fi
        done < submodules_config.txt
        git submodule init
        git submodule update --remote --merge
      else
        echo "[INFO] 未找到 submodules_config.txt"
      fi
      ;;
    3)
      echo "[INFO] 更新主專案..."
      if ! git remote get-url upstream > /dev/null 2>&1; then
        read -p "請輸入 upstream URL: " upstreamURL
        git remote add upstream "$upstreamURL"
      fi
      git fetch upstream
      git pull upstream main --allow-unrelated-histories

      if [ -f "submodules_config.txt" ]; then
        while read -r line; do
          set -- $line
          subPath="$1"
          originURL="$2"
          upstreamURL="$3"
          if [ -d "$subPath" ]; then
            cd "$subPath"
            git remote set-url origin "$originURL"
            if ! git remote get-url upstream > /dev/null 2>&1; then
              git remote add upstream "$upstreamURL"
            else
              git remote set-url upstream "$upstreamURL"
            fi
            git checkout main
            git fetch upstream
            git pull upstream main
            cd ..
          fi
        done < submodules_config.txt
      fi
      ;;
    4)
      read -p "請輸入 commit 訊息（預設：更新）: " commitMsg
      [ -z "$commitMsg" ] && commitMsg="更新"
      timestamp=$(date "+%Y-%m-%d_%H-%M")
      echo "Commit Log - $timestamp" > commit_log.txt
      echo "-------------------------" >> commit_log.txt

      if [ -f "submodules_config.txt" ]; then
        echo "[INFO] 推送 submodules..."
        while read -r line; do
          set -- $line
          subPath="$1"
          if [ -d "$subPath" ]; then
            cd "$subPath"
            git add .
            git commit -m "$commitMsg - $timestamp" 2>/dev/null
            git push origin main
            echo "[submodule] $subPath 已提交：$commitMsg - $timestamp" >> ../commit_log.txt
            cd ..
          fi
        done < submodules_config.txt
      fi

      echo "[INFO] 提交主專案"
      git add .
      git commit -m "$commitMsg - $timestamp" 2>/dev/null
      git push origin main
      echo "[main] 主專案已提交：$commitMsg - $timestamp" >> commit_log.txt
      ;;
    X|x)
      echo "👋 Bye!"
      break
      ;;
    *)
      echo "❌ 無效的選項，請重新選擇"
      ;;
  esac
done