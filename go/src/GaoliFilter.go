// GetSeq
package main

import (
	"bufio"
	"bytes"
	"flag"
	//	"fmt"
	"lpp"
	"os"
)

func main() {

	/* Get Option Parser
	 */
	blast := flag.String("b", "", "Read1")

	output := flag.String("o", "", "Output")

	flag.Parse()

	//	fmt.Println(data1, data2)
	/* Generate Output
	 */
	OUTPUT, err := os.Create(*output + ".stats")

	if err != nil {
		panic("Can not Create Result File!!")
	}
	BufOUTPUT := bufio.NewWriterSize(OUTPUT, 9999)
	dataIO := lpp.GetBlockRead( *blast, "\n", false, 10000000)
	for {

		line, err := dataIO.Next()
		if err != nil {
			break
		}

		for _, e_sub := range [][]byte{[]byte("Bacteria"), []byte("Chordata"), []byte("Arthropoda"), []byte("Viruses"), []byte("Streptophyta")} {
			if bytes.Contains(line, e_sub) {
				BufOUTPUT.Write(line)
				break
			}
		}

	}

	BufOUTPUT.Flush()

}
