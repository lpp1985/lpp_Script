package main

import (
	"bytes"
	"flag"
	"lpp"
)

func main() {

	inputdata := flag.String("i", "", "Input")
	start := flag.Int("s", 0, "Start Location")

	end := flag.Int("e", 150, "Start Location")
	outdata := flag.String("o", "output.fastq", "Start Location")

	end := flag.Int("e", 150, "End Location")
	outdata := flag.String("o", "output.fastq", "Output")

	flag.Parse()
	OUTPUT, _ := lpp.GetOuput(*outdata, 1000000)
	FQ := lpp.Fastq{File: *inputdata}
	for {
		a, b, c, d, err := FQ.Next()
		seq := bytes.TrimSpace(b)[*start:*end]
		qual := bytes.TrimSpace(d)[*start:*end]
		OUTPUT.Write(a)
		OUTPUT.Write(seq)
		OUTPUT.WriteString("\n")
		OUTPUT.Write(c)
		OUTPUT.Write(qual)
		OUTPUT.WriteString("\n")
		if err != nil {
			break
		}

	}

}
