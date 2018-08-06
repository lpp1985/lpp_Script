package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"lpp"
	"os"
)

func main() {
	err := recover()
	if err != nil {
		fmt.Println(err)
	}
	input := flag.String("i", "", "input")

	output := flag.String("o", "", "Output")
	list_number := flag.Int("n", 1, "list Number")

	flag.Parse()

	OUTPUTHANDLE, err := os.Create(*output)
	defer OUTPUTHANDLE.Close()
	OUTPUTBUF := bufio.NewWriterSize(OUTPUTHANDLE, 10000)
	//	OUTPUTBUF.WriteString("Name\tLength\n")
	defer OUTPUTBUF.Flush()
	if err != nil {
		panic("Output not Exist!!")
	}

	INPUTIO := lpp.GetBlockRead(*input, "\n", false, 10000000)
	all_has := make(map[string]string)
	for {
		line, err := INPUTIO.Next()

		name := string(bytes.Fields(line)[*list_number-1])
		_, ok := all_has[name]
		if !ok {
			OUTPUTBUF.Write(line)
		}
		all_has[name] = ""

		if err != nil {
			break
		}
	}

}
