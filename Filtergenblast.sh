grep "Gene:ID=" $1| sed -r "s/Gene:ID=//g"| sed -r "s/\|/\t/g"| sed -r "s/_All\S+//g" 
