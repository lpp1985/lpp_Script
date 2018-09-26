// KmerFinder
package main

import (
	"bytes"
	"flag"
	"fmt"
	"lpp"
)

var kmer_number map[string]int = make(map[string]int)
var kmer_graph map[string]map[string]int = make(map[string]map[string]int)
var total_mer int
var kmer_spices int

func Kmer_Count(file_handle string, kmer int) {
	F := lpp.Fastq{File: file_handle}
	for {

		_, seq, _, _, err := F.Next()
		seq = bytes.TrimSpace(seq)
		for j := 0; j <= len(seq)-kmer-1; j++ {
			i := j + 1
			mer1 := string(seq[j : j+kmer])
			mer2 := string(seq[i : i+kmer])

			_, ok1 := kmer_number[mer1]
			if !ok1 {
				kmer_spices += 1
				kmer_number[mer1] = 1
			} else {
				kmer_number[mer1] += 1
			}
			_, ok2 := kmer_number[mer2]
			if !ok2 {
				kmer_spices += 1
				kmer_number[mer2] = 1
			} else {
				kmer_number[mer2] += 1
			}
			_, ok_graph := kmer_graph[mer1]
			if !ok_graph {
				kmer_graph[mer1] = make(map[string]int)
				kmer_graph[mer1][mer2] = 1
			} else {
				kmer_graph[mer1][mer2] += 1
			}
		}
		if err != nil {
			break
		}
	}
}

func main() {
	read1 := flag.String("1", "", "read1")
	read2 := flag.String("2", "", "read2")
	output := flag.String("o", "", "Output")
	kmer := flag.Int("k", 41, "kmer")
	flag.Parse()
	Output, _ := lpp.GetOuput(*output, 10000)
	fmt.Println("read1 Start!!")
	Kmer_Count(*read1, *kmer)
	fmt.Println("read2 Start!!")
	Kmer_Count(*read2, *kmer)
	for m1, h2 := range kmer_graph {
		for m2, number := range h2 {
			Output.WriteString(m1 + "\t" + m2 + "\t" + fmt.Sprintf("%d", number) + "\n")
		}
	}
}
