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

	/* Get Option Parser
	 */
	input := flag.String("i", "", "Input Data from ftp://ftp.pir.georgetown.edu/databases/idmapping/idmapping.tb.gz ")
	id_uniprot := flag.String("u", "", "id-->uniprotID Mapping File")
	uniprot_go := flag.String("g", "", "uniprotID -->goid Mapping File")
	gi := flag.String("k", "", "UniprotID-->GI Mapping File")

	flag.Parse()
	/* Generate Output
	 */
	Uniprot_File, err := os.Create(*id_uniprot)
	if err != nil {
		panic("Can not Create id-->uniprotID  File!!")
	}
	defer Uniprot_File.Close()
	BufUniprot := bufio.NewWriterSize(Uniprot_File, 9999999)

	GO_File, err := os.Create(*uniprot_go)
	if err != nil {
		panic("Can not Create uniprotID -->goid File!!")
	}
	defer GO_File.Close()
	BufGO := bufio.NewWriterSize(GO_File, 9999999)

	GI_File, err := os.Create(*gi)
	if err != nil {
		panic("Can not Create UniprotID-->GI Mapping File!!")
	}

	defer GI_File.Close()
	BufGI := bufio.NewWriterSize(GI_File, 9999999)

	/*
		open idmapping.tb Input

	*/
	dataIO := lpp.GetBlockRead(*input, "\n", false, 100000000)
	uniprot_id := 0
	for {

		line, err := dataIO.Next()
		if err != nil {
			break
		}
		uniprot_id += 1
		line_l := bytes.Split(line, []byte("\t"))
		BufUniprot.WriteString(fmt.Sprintf("%d\t%s\n", uniprot_id, line_l[0]))
		if len(line_l[7]) != 0 {
			for _, e_go := range bytes.Split(line_l[7], []byte("; ")) {

				BufGO.WriteString(fmt.Sprintf("%d\t%s\n", uniprot_id, e_go))

			}
		}
		for _, e_gi := range bytes.Split(line_l[4], []byte("; ")) {
			if e_gi != nil {
				BufGI.WriteString(fmt.Sprintf("%d\t%s\n", uniprot_id, e_gi))
			}

		}

	}
	defer BufGO.Flush()
	defer BufGI.Flush()
	defer BufUniprot.Flush()

}
