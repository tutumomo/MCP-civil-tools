#!/bin/bash

echo "========================================"
echo "[Menu] Git 操作選擇"
echo "========================================"
echo "[0] Initialize remote & submodules (already cloned)"
echo "[1] Initialize submodules (with check)"
echo "[2] Pull updates from upstream (main + submodules)"
echo "[3] Push local changes to GitHub (main + submodules)"
echo "[X] Exit"
echo "========================================"

read -p "Enter your choice: " choice

function submodule_init_line() {
    subPath="$1"
    originURL="$2"
    upstreamURL="$3"

    existing_url=$(git config -f .gitmodules --get submodule."$subPath".url)
    if [ -z "$existing_url" ]; then
        echo "[INFO] Adding submodule: $subPath"
        git submodule add "$originURL" "$subPath"
    else
        echo "[INFO] Submodule $subPath already exists. Skipping."
    fi
}

case "$choice" in
  0)
    echo "[INFO] Checking upstream remote..."
    git remote get-url upstream >/dev/null 2>&1
    if [ $? -ne 0 ]; then
        read -p "Enter upstream URL for main repo: " upstreamURL
        git remote add upstream "$upstreamURL"
    else
        echo "[INFO] Upstream already exists. Skipping."
    fi

    if [ -f "submodules_config.txt" ]; then
        echo "[INFO] Initializing submodules..."
        while read -r line; do
            set -- $line
            submodule_init_line "$1" "$2" "$3"
        done < submodules_config.txt
        git submodule init
        git submodule update --remote --merge
    else
        echo "[INFO] No submodules_config.txt found. Skipping."
    fi
    ;;

  1)
    if [ -f ".gitmodules" ]; then
        read -p "[INFO] Submodules already exist. Re-initialize? (Y/N): " redoInit
        if [[ ! "$redoInit" =~ ^[Yy]$ ]]; then
            echo "[INFO] Submodule init canceled."
            exit 0
        fi
    fi

    if [ -f "submodules_config.txt" ]; then
        echo "[INFO] Initializing submodules from submodules_config.txt..."
        while read -r line; do
            set -- $line
            submodule_init_line "$1" "$2" "$3"
        done < submodules_config.txt
        git submodule init
        git submodule update --remote --merge
    else
        echo "[INFO] No submodules_config.txt found. Skipping."
    fi
    ;;

  2)
    echo "[INFO] Pulling from upstream..."
    git remote get-url upstream >/dev/null 2>&1
    if [ $? -ne 0 ]; then
        read -p "Enter upstream URL for main repo: " upstreamURL
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
                git remote get-url upstream >/dev/null 2>&1
                if [ $? -ne 0 ]; then
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

  3)
    read -p "Enter commit message (default: 更新): " commitMsg
    [ -z "$commitMsg" ] && commitMsg="更新"
    timestamp=$(date "+%Y-%m-%d_%H-%M")
    echo "Commit Log - $timestamp" > commit_log.txt
    echo "--------------------------" >> commit_log.txt

    if [ -f "submodules_config.txt" ]; then
        echo "[INFO] Pushing submodules..."
        while read -r line; do
            set -- $line
            subPath="$1"
            if [ -d "$subPath" ]; then
                cd "$subPath"
                git add .
                git commit -m "$commitMsg - $timestamp" 2>/dev/null
                git push origin main
                echo "[submodule] $subPath committed: $commitMsg - $timestamp" >> ../commit_log.txt
                cd ..
            fi
        done < submodules_config.txt
    fi

    echo "[INFO] Committing main repo"
    git add .
    git commit -m "$commitMsg - $timestamp" 2>/dev/null
    git push origin main
    echo "[main] Main repo committed: $commitMsg - $timestamp" >> commit_log.txt
    ;;

  X|x)
    echo "👋 Exit"
    exit 0
    ;;

  *)
    echo "❌ Invalid choice"
    ;;
esac