
CLS
DEFDBL A-Z
PI = 4 * ATN(1)
PP = 180/PI
PR = PI/180
JFS = 1400*0.875
PP$ = ""
I% = 0
PRINT
FF$ = ""
INPUT "結 果 存 放 檔 名 ";FF$
PRINT
IF FF$ <> "" THEN
   METHOD% = 0
   ON ERROR GOTO ERRHANDLE
   INPUT "原 有 檔 案 是 否 刪 除 (Y:刪除  N:保留)";ANS$
   PRINT
   IF ANS$ = "Y" THEN
      IF METHOD% = 0 THEN OPEN "O",#4,FF$
   ELSE
      IF METHOD% = 0 THEN OPEN "A",#4,FF$
   END IF
   ON ERROR GOTO 0
END IF
ON ERROR GOTO ERRHANDLE
PP$ = ""
GOTPP :
INPUT "是 否 列 印 計 算 數 據 (Y/N)";PP$
IF PP$="y" THEN PP$ = "Y"
INPUT "  工程名稱 : ";NA$
PRINT
INPUT "  土    重 ( ｒ ) : ";R
PRINT
IF PP$ = "Y" THEN
   LPRINT
   LPRINT "  工程名稱 : ";NA$;TAB(64);"日期: ";DATE$
   LPRINT "  土  重 (ｒ) = ";
   LPRINT USING "#.###";R
   LPRINT "==============================================================================="
   LPRINT "      高 度 溝壁傾角 土方傾角 安息角 土壓係數 土壓力  彎  距  有效厚度 鋼 筋 量"
   LPRINT "編號   H(m)   ｍ(m)    (ｉ)    (ψ)   (Ｋａ)   (Ｐ)    (Ｍ)    ｄ(m)   Ａｓ()"
   LPRINT "---- ------ -------- -------- ------ -------- ------- ------- -------- --------"
END IF
IF FF$ <> "" THEN
   PRINT #4,""
   PRINT #4,"  工程名稱 : ";NA$;TAB(64);"日期: ";DATE$
   PRINT #4,"  土  重 (ｒ) = ";
   PRINT #4,USING "#.###";R
   PRINT #4,"==============================================================================="
   PRINT #4,"      高 度 溝壁傾角 土方傾角 安息角 土壓係數 土壓力  彎  距  有效厚度 鋼 筋 量"
   PRINT #4,"編號   H(m)   ｍ(m)    (ｉ)    (ψ)   (Ｋａ)   (Ｐ)    (Ｍ)    ｄ(m)   Ａｓ()"
   PRINT #4,"---- ------ -------- -------- ------ -------- ------- ------- -------- --------"
END IF
AGAIN :
INPUT "  編    號 ( No ) : ";PN$
PRINT
INPUT "  高    度 ( Ｈ ) : ";H
PRINT
INPUT "  溝壁傾角 ( ｍ ) : ";J
PRINT
INPUT "  土方傾角 ( ｉ ) : ";I
PRINT
INPUT "  安 息 角 ( ψ ) : ";Y
PRINT
INPUT "  有效厚度 ( ｄ ) : ";D
PRINT
IF PP$ = "Y" THEN
   LPRINT TAB(1);PN$;
   LPRINT TAB(7);
   LPRINT USING "#.###";H;
   LPRINT TAB(15);
   LPRINT USING "#.###";J;
   LPRINT TAB(24);
   LPRINT USING "##.##";I;
   LPRINT TAB(32);
   LPRINT USING "##.##";Y;
END IF
IF FF$ <> "" THEN
   PRINT #4,TAB(1);PN$;
   PRINT #4,TAB(7);
   PRINT #4,USING "#.###";H;
   PRINT #4,TAB(15);
   PRINT #4,USING "#.###";J;
   PRINT #4,TAB(24);
   PRINT #4,USING "##.##";I;
   PRINT #4,TAB(32);
   PRINT #4,USING "##.##";Y;
END IF
PRINT "請 稍 待 ......"
PRINT
JA = ATN(J)
IA = I*PR
YA = Y*PR
IF (SIN(YA)*SIN(YA-IA)) / (COS(JA+IA)*COS(J)) < 0 THEN
   PRINT " 根 據 公 式 此 資 料 不 能 計 算 "
   GOTO ISCONT
END IF
KA = (COS(YA+JA))^2 / ( COS(JA)^2 * ( 1+SQR( (SIN(YA)*SIN(YA-IA)) / (COS(JA+IA)*COS(J)) ) )^2)
MP = R*H^2*KA/(2*COS(JA))
MB = R*H^3*KA/(6*COS(JA))
PRINT "  土壓係數 (Ｋａ) = ";
PRINT USING "###.#######";KA
PRINT
PRINT "  土 壓 力 ( Ｐ ) = ";
PRINT USING "###.#######";MP
PRINT
PRINT "  彎    距 ( Ｍ ) = ";
PRINT USING "###.###";MB
PRINT
SA = CLNG(MB / (JFS*D) * 1000000)/1000
PRINT "  鋼筋需求量 (Ａｓ) = ";
PRINT USING "###.###";SA
IF PP$ = "Y" THEN
   LPRINT TAB(39);
   LPRINT USING "#.####";KA;
   LPRINT TAB(47);
   LPRINT USING "###.###";MP;
   LPRINT TAB(55);
   LPRINT USING "###.###";MB;
   LPRINT TAB(64);
   LPRINT USING "##.###";D;
   LPRINT TAB(72);
   LPRINT USING "###.###";SA
END IF
IF FF$ <> "" THEN
   PRINT #4,TAB(39);
   PRINT #4,USING "#.####";KA;
   PRINT #4,TAB(47);
   PRINT #4,USING "###.###";MP;
   PRINT #4,TAB(55);
   PRINT #4,USING "###.###";MB;
   PRINT #4,TAB(64);
   PRINT #4,USING "##.###";D;
   PRINT #4,TAB(72);
   PRINT #4,USING "###.###";SA
END IF
ISCONT :
PRINT
INPUT "是 否 繼 續 ? (Y/N)";ANS$
PRINT
IF ANS$ = "Y" OR ANS$ = "y" THEN AGAIN
CLOSE
END

ERRHANDLE :
KK% = 0
SELECT CASE ERR
   CASE 2
        PRINT "資 料 輸 入 錯 誤 , 請 重 新 輸 入 ...."
        KK% = 1
   CASE 11
        PRINT "計 算 有 誤 , 請 檢 查 資 料"
        KK% = 2
   CASE 25, 27
        PRINT "列表機未裝置妥當 , 裝好後按任一鍵繼續 或 按[Ａ]鍵結束"
        WHILE NOT INSTAT
        WEND
        S$ = INKEY$
        IF S$ = "A" OR S$="a" THEN CLOSE: END
        KK% = 1
   CASE 64, 53
        PRINT "檔 名 輸 入 錯 誤 , 請 重 新 輸 入"
        INPUT "結 果 存 放 檔 名 ";FF$
        PRINT
        IF FF$ = "" THEN
           METHOD% = 1
        END IF
        KK% = 1
END SELECT
IF KK% = 0 THEN
   PRINT ERR
   PRINT "程 式 無 法 處 理"
   PRINT
   CLOSE
   END
END IF
IF KK% = 2 THEN
   IF PP$="Y" THEN
      LPRINT
   END IF
   RESUME AGAIN
END IF
RESUME 0

