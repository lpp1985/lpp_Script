git add . -A 
git commit -m 123 
git pull -f Roomnas master:$1
git submodule init 
git submodule update
