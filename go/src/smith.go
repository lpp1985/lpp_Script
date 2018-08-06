package main

import (
	"bufio"
	"bytes"
	"flag"
	//	"fmt"
	"lpp"
	"os"
)

func Select(seq string) string {
	var max float64
	result := ""

	for each_index, name := range index_seq {
		score := lpp.SmithWaterman(seq, each_index)
		if score > max {
			max = score
			result = name
		}
	}
	return result

}

var index_seq map[string]string = make(map[string]string)

type PairFastq struct {
	Pair1 *lpp.IO
	Pair2 *lpp.IO
}

func (Self *PairFastq) Next() ([4][]byte, [4][]byte, error) {
	data1 := [4][]byte{}
	data2 := [4][]byte{}
	i := 0
	for {
		i += 1
		line1, err1 := Self.Pair1.Next()
		if err1 != nil {
			return data1, data2, err1
		}
		data1[i-1] = line1
		line2, _ := Self.Pair2.Next()
		data2[i-1] = line2
		if i == 4 {
			return data1, data2, err1
		}

	}

}

func main() {

	index := flag.String("i", "", "Index")
	outputdir := flag.String("o", "", "Output")
	flag.Parse()
	fastq1 := flag.String("1", "", "fastq1")
	fastq2 := flag.String("2", "", "fastq1")

	FQ1IO := lpp.GetBlockRead(*fastq1, "\n", false, 10000000)

	FQ2IO := lpp.GetBlockRead(*fastq2, "\n", false, 10000000)
	FQ1Q2IO := PairFastq{}
	FQ1Q2IO.Pair1 = &FQ1IO
	FQ1Q2IO.Pair2 = &FQ2IO
	if _, err := os.Stat(*outputdir); os.IsNotExist(err) {
		os.Mkdir(*outputdir, os.ModePerm)
	}
	*outputdir = *outputdir + "/"

	DATAIO := lpp.GetBlockRead(*index, "\n", false, 10000000)

	seq_result := make(map[string][2]*bufio.Writer)
	for {
		line, err := DATAIO.Next()
		if err != nil {
			break
		}
		line = bytes.TrimSuffix(line, []byte("\n"))
		line = bytes.TrimPrefix(line, []byte("\n"))
		name := string(bytes.Fields(line)[0])
		seq := string(bytes.Fields(bytes.SplitN(line, []byte("\t"), 2)[1])[0])

		index_seq[seq] = name
		resultFile1, err1 := os.Create(*outputdir + name + ".R1.fastq")
		if err1 != nil {
			panic("Fastq not Exist " + *outputdir + name + ".R1.fastq")

		}
		OUTPUTBUF1 := bufio.NewWriterSize(resultFile1, 10000)
		resultFile12, err2 := os.Create(*outputdir + name + ".R2.fastq")
		if err2 != nil {
			panic("Fastq not Exist " + *outputdir + name + ".R2.fastq")

		}
		OUTPUTBUF2 := bufio.NewWriterSize(resultFile12, 10000)
		seq_result[string(name)] = [2]*bufio.Writer{OUTPUTBUF1, OUTPUTBUF2}

	}
	for {
		fq1data, fq2data, err := FQ1Q2IO.Next()
		if err != nil {
			break
		}
		seq_f2 := string(fq2data[1])
		indexSEQ := seq_f2[23 : 23+6]
		result_index, ok := index_seq[indexSEQ]
		var result_handle [2]*bufio.Writer
		if ok {
			result_handle = seq_result[result_index]

		} else {
			final_index := Select(indexSEQ)
			result_handle = seq_result[final_index]
		}

		for _, cont := range fq1data {
			result_handle[0].Write(cont)
		}
		for _, cont := range fq2data {
			result_handle[1].Write(cont)
		}

	}

}
