package main

import (
	//	"bytes"
	"fmt"
	"lpp"
	"os"
)

func main() {
	FASTQ := lpp.Fastq{File: os.Args[1]}

	//	if FASTA.IO.BlockTag == nil {
	//		fmt.Println("Blank!!")
	//	}

	for {
		name, seq, name2, qual, err := FASTQ.Next()
		fmt.Println(string(name[len(name)-1]))
		fmt.Println(string(seq[len(seq)-1]))

		fmt.Println(string(name2[len(name2)-1]))
		fmt.Println(string(qual[len(qual)-1]))
		if err != nil {
			break
		}

	}
}
