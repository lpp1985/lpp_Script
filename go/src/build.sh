for i in *.go ; do go build  -o ../bin/${i%.go} $i ; done  
