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
	fasta := flag.String("f", "", "fasta")

	output := flag.String("o", "", "Output")

	flag.Parse()

	OUTPUTHANDLE, err := os.Create(*output)
	defer OUTPUTHANDLE.Close()
	OUTPUTBUF := bufio.NewWriterSize(OUTPUTHANDLE, 10000)
	OUTPUTBUF.WriteString("Name\tLength\n")
	defer OUTPUTBUF.Flush()
	if err != nil {
		panic("Output not Exist!!")
	}

	FASTAIO := lpp.GetBlockRead(*fasta, "\n>", false, 10000000)
	for {
		line, err := FASTAIO.Next()

		line = bytes.TrimSuffix(line, []byte(">"))
		line = bytes.TrimPrefix(line, []byte(">"))
		name := bytes.Fields(line)[0]
		seq := bytes.SplitN(line, []byte("\n"), 2)[1]
		seq = bytes.Replace(seq, []byte("\n"), []byte(""), -1)
		length := len(seq)

		output_byte := name
		output_byte = append(output_byte, []byte("\t")...)
		output_byte = append(output_byte, []byte(string(length))...)
		output_byte = append(output_byte, []byte("\n")...)
		OUTPUTBUF.Write(output_byte)
		if err != nil {
			break
		}
	}

}
