chmod -R 755 ./*.*
#dos2unix *.*
#GIT_TRACE=1
git add . -A
git commit -m 123
#git push -f Homenas pc:master
git push -f Roomnas $1:master&
git push -f USA $1:master&
git push -f GIT $1:master&
git push -f HG $1:master&
