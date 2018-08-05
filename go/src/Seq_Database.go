// N50
package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	. "lpp"
	"os"
	//	"sort"
)

func main() {

	//	var length_Ddict map[int]map[string]string = make(map[int]map[string]string)

	file := flag.String("i", "", "input Fasta!")
	output := flag.String("o", "", "Output!")
	flag.Parse()
	if *file == "" {
		os.Exit(1)
	}
	fasta := new(Block_Reading)
	rawfile,_:=os.Open(*file)
	fasta.File = rawfile
	fasta.Blocktag = "\n>"
	fasta_handle:= fasta.Read()
	WRHANDLE, _ := os.Create(*output)
	RESULT := bufio.NewWriterSize(WRHANDLE,1000000)

	for {
		line, err := fasta_handle.Next()
		data := bytes.SplitN(line, []byte("\n"), 2)
		seq := data[1]
		title := data[0]
		cache := bytes.Fields(title)
		name := cache[0]
		//fmt.Println(string(name))
		//os.Exit(1)
		annotation := string(title)
		if name[0] == '>' {
			name = name[1:]
		}

		//		name := string(data[0])
		seq = bytes.Replace(seq, []byte("\n"), []byte(""), -1)
		seq_length := len(seq)
		
		RESULT.WriteString(fmt.Sprintf("%s\t%s\t%d\n", name, annotation, seq_length))

		//		_, ok := length_Ddict[all_length][name]
		//		if !ok {
		//			length_Ddict[all_length] = make(map[string]string)
		//			length_Ddict[all_length][name] = ""
		//		}
		if err != nil {
			break
		}

	}
	RESULT.Flush()
}
