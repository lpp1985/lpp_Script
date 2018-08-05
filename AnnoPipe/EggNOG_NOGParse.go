// GetSeq
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
	var all_has_hash map[string]string = make(map[string]string)

	/* Get Option Parser
	 */
	input := flag.String("i", "", "Input Annotation Data")
	annotation := flag.String("a", "", "COG annotation")
	category := flag.String("c", "", "COG Category")
	raw := flag.String("r", "", "Input Member Data")
	mem := flag.String("m", "", "NOG Member Result")
	fasta := flag.String("f", "", "Fasta Input")
	sequence := flag.String("s", "", "Fasta output")
	flag.Parse()
	/* Generate Output
	 */
	Mem_File, err := os.Create(*mem)
	if err != nil {
		panic("Can not Create Annotation File!!")
	}
	defer Mem_File.Close()
	BufMem := bufio.NewWriterSize(Mem_File, 9999)

	Annotation_File, err := os.Create(*annotation)
	if err != nil {
		panic("Can not Create Annotation File!!")
	}
	defer Annotation_File.Close()
	BufAnnotation := bufio.NewWriterSize(Annotation_File, 9999)

	Category_File, err := os.Create(*category)
	if err != nil {
		panic("Can not Create Category File!!")
	}

	defer Category_File.Close()
	BufCategory := bufio.NewWriterSize(Category_File, 9999)

	/*
		open Fasta Input

	*/
	dataIO, err := lpp.GetBlockRead(*input, "\n", false, 100000000)
	if err != nil {
		panic("Input Input Error")
	}

	for {

		line, err := dataIO.Next()
		if err != nil {
			break
		}
		line_l := bytes.Split(line, []byte("\t"))
		BufAnnotation.WriteString(fmt.Sprintf("%s\t%s", line_l[1], line_l[5]))
		if err != nil {
			break
		}
		for _, e_cat := range bytes.Split(line_l[4], []byte("")) {
			BufCategory.WriteString(fmt.Sprintf("%s\t%s\n", line_l[1], e_cat))
		}
	}
	defer BufCategory.Flush()
	defer BufAnnotation.Flush()
	memIO, err := lpp.GetBlockRead(*raw, "\n", false, 100000000)
	for {

		line, err := memIO.Next()
		if err != nil {
			break
		}
		line_l := bytes.Split(line, []byte("\t"))
		all_mem := bytes.Replace(line_l[5], []byte("\n"), []byte(""), -2)

		for _, e_mem := range bytes.Split(all_mem, []byte(",")) {
			BufMem.WriteString(fmt.Sprintf("%s\t%s\n", e_mem, line_l[1]))
			all_has_hash[string(e_mem)] = ""
		}
	}
	defer BufMem.Flush()
	Seq_File, err := os.Create(*sequence)
	if err != nil {
		panic("Can not Create Sequence File!!")
	}
	defer Seq_File.Close()
	BufSeq := bufio.NewWriterSize(Seq_File, 9999)

	fastaIO, err := lpp.GetBlockRead(*fasta, "\n>", false, 1000000000)
	for {
		line, err := fastaIO.Next()
		if err != nil {
			break
		}
		line = bytes.TrimSuffix(line, []byte(">"))
		name := bytes.SplitN(line, []byte("\n"), 2)[0]
		name = bytes.Fields(name)[0]
		_, has := all_has_hash[string(name)]
		if has {
			//			fmt.Println(i)
			line = append([]byte(">"), line...)
			BufSeq.Write(line)

		}
	}
	defer BufSeq.Flush()

}
