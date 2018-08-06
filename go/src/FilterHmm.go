package main

import (
	"bufio"
	//	"bytes"
	"flag"
	"fmt"
	"lpp"
	"os"
	"regexp"
	"strings"
)

func main() {
	err := recover()
	if err != nil {
		fmt.Println(err)
	}
	hmmer := flag.String("f", "", "hmmer")
	list := flag.String("l", "", "list")
	number := flag.Int("n", 2, "number")
	output := flag.String("o", "", "Output")

	flag.Parse()

	OUTPUTHANDLE, err := os.Create(*output)
	defer OUTPUTHANDLE.Close()
	OUTPUTBUF := bufio.NewWriterSize(OUTPUTHANDLE, 10000)
	//	OUTPUTBUF.WriteString("Name\tLength\n")
	defer OUTPUTBUF.Flush()
	if err != nil {
		panic("Output not Exist!!")
	}
	HMMERHANDLE, err := os.Open(*hmmer)
	LISTHANDLE, err := os.Open(*list)
	if err != nil {
		panic("Database error")
	}
	defer HMMERHANDLE.Close()
	if err != nil {
		panic("Fasta not Exist")
	}
	HMMERIO := lpp.GetBlockRead(HMMERHANDLE, "\n//\n", false, 10000000)
	re := regexp.MustCompile(`\nNAME\s+([^\.]+)`)

	dataIO := lpp.GetBlockRead(LISTHANDLE, "\n", false, 1000000)

	/* need_hash Generate
	 */
	raw_hash := new(lpp.File_dict)
	raw_hash.File_IO = dataIO
	raw_hash.Header = false
	need_hash := raw_hash.Read(*number, *number)
	fmt.Println(len(need_hash))
	already := make(map[string]string)
	for {

		line, err := HMMERIO.Next()

		if err != nil {
			break
		}
		name := string(re.FindSubmatch(line)[1])
		_, has := need_hash[name]
		if has {
			already[name] = ""
			OUTPUTHANDLE.WriteString(string(line))
		}

	}
	for e_id, _ := range need_hash {
		_, kk := already[e_id]
		if !kk {
			fmt.Println(e_id)

		}
	}

}
