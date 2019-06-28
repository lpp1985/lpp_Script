// KmerFinder
package main

import (
	"bytes"
	"flag"
	"fmt"
	"lpp"
)

func Int_Byte(data int64, kmer int64) []byte {
	var result = make([]byte, kmer)
	for i := int64(0); i <= kmer-1; i++ {
		result[i] = 'A'
	}
	for i := int64(1); i <= kmer-1; i++ {

		yu := data % 4
		data = data / 4
		result[kmer-i] = complement[uint8(yu)]
		if data == 0 {
			break
		}

	}
	return result
}

var complement = [256]uint8{
	0: 'A', 1: 'T',
	2: 'C', 3: 'G',
	'A': 0, 'a': 0,
	'C': 2, 'c': 2,
	'G': 3, 'g': 3,
	'T': 1, 't': 1,
}

func exponent(a, n int64) int64 {
	result := int64(1)
	for i := n; i > 0; i >>= 1 {
		if i&1 != 0 {
			result *= a
		}
		a *= a
	}
	return result
}
func Byte_Int(str []byte) int64 {
	l := int64(len(str))
	var num int64
	for i := l - 1; i >= 0; i-- {
		data := str[i]
		loc := int64(complement[data])
		if l-i-1 == 0 {
			num += loc
		} else {
			num += exponent(4, (l-i-1)) * loc
		}
	}
	return num

}

var kmer_number map[string]int = make(map[string]int)
var kmer_graph map[string]map[string]int = make(map[string]map[string]int)

func Kmer_Count(file_handle string) {
	F := lpp.Fasta{File: file_handle}
	for {

		_, seq, err := F.Next()
		seq = bytes.TrimSpace(seq)
		is_n := bytes.Contains(seq, []byte("N"))
		if is_n {
			continue
		}

		kmer_number[string(seq)] = 1

		if err != nil {
			break
		}
	}

}

func main() {
	read1 := flag.String("f", "", "read1")
	kmer := flag.Int("k", 41, "kmer")
	output := flag.String("o", "", "Output")

	flag.Parse()
	Output, _ := lpp.GetOuput(*output, 10000)
	fmt.Println("read1 Start!!")
	Kmer_Count(*read1)

	fmt.Println(len(kmer_number))
	for raw_mer, _ := range kmer_number {

		mer_suff := raw_mer[1:*kmer]
		for _, i := range [4]string{"A", "T", "C", "G"} {
			next_mer := mer_suff + i
			_, ok := kmer_number[next_mer]
			if ok {
				_, ok2 := kmer_graph[raw_mer]
				if !ok2 {
					kmer_graph[raw_mer] = make(map[string]int)
				}
				kmer_graph[raw_mer][next_mer] = 0
			}
		}
	}

	for mer1, h1 := range kmer_graph {
		for mer2, number := range h1 {
			Output.WriteString(fmt.Sprintf("%s\t%s\t%d\n", mer1, mer2, number))
		}
	}

}
