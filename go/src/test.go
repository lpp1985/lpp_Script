package main

import (
	"bytes"
	"flag"
	"fmt"

	//"fmt"
	"lpp"
)

func main() {

	input := flag.String("i", "./", "Input File")

	flag.Parse()
	raw_file := lpp.GetBlockRead(*input, "\n", false, 10000)
	var All_list []string
	for {
		line, err := raw_file.Next()
		line = bytes.TrimSpace(line)
		All_list = append(All_list, string(line))
		fmt.Println(string(line))
		if err != nil {
			break
		}

	}
	RAW := lpp.Fasta{File: *input}
	for {
		title, _, err := RAW.Next()

		fmt.Println(string(title))
		if err != nil {
			break
		}
	}
}
