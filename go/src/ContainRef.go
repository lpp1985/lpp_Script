package main

import (
	"bytes"
	"flag"
	//	"fmt"
	"lpp"
)

func main() {
	input := flag.String("i", "./", "Input Fasta")
	database := flag.String("l", "./", "TransGene Sequence")
	output := flag.String("o", "./", "Output")
	flag.Parse()
	OUTPUT, _ := lpp.GetOuput(*output, 1000)
	RAWIO := lpp.Fasta{File: *input}
	DATAIO := lpp.Fasta{File: *database}
	all_need := [][]byte{}
	for {
		_, s, err := DATAIO.Next()
		s = bytes.Fields(bytes.ToUpper(s))[0]
		s1 := lpp.RevComplement(s)
		all_need = append(all_need, s1[len(s1)-41:])
		all_need = append(all_need, s[len(s)-41:])
		if err != nil {
			break
		}
	}
	//	for _, data := range all_need {
	//		fmt.Println(string(data))
	//	}

	for {
		t, s, err := RAWIO.Next()
		s = bytes.ToUpper(s)
		s1 := lpp.RevComplement(s)
		t = bytes.Fields(t)[0][1:]
		if bytes.Contains(s, all_need[0]) {
			OUTPUT.Write(t)
			OUTPUT.WriteString("+\n")
			break
		} else {
			if bytes.Contains(s1, all_need[0]) {
				OUTPUT.Write(t)
				OUTPUT.WriteString("-\n")
				break
			}
		}

		if err != nil {
			break
		}
	}
	//	OUTPUT.WriteString("Head Found!!!\n")
	RAWIO = lpp.Fasta{File: *input}
	for {
		t, s, err := RAWIO.Next()
		s = bytes.ToUpper(s)
		s1 := lpp.RevComplement(s)
		t = bytes.Fields(t)[0][1:]
		if bytes.Contains(s, all_need[1]) {
			OUTPUT.Write(t)
			OUTPUT.WriteString("+\n")
			break
		} else {
			if bytes.Contains(s1, all_need[1]) {
				OUTPUT.Write(t)
				OUTPUT.WriteString("-\n")
				break
			}
		}

		if err != nil {
			break
		}
	}

}
