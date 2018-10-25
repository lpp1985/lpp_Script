// KmerFinder
package main

import (
	"bytes"
	"flag"
	"fmt"
	"lpp"
	"math/big"
	"sort"
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

var kmer_number map[int64]int = make(map[int64]int)
var kmer_graph map[string]map[string]int = make(map[string]map[string]int)
var total_mer int
var kmer_spices int
var repre_mer map[string]string = make(map[string]string)

func Kmer_Count(file_handle string, kmer int) {
	F := lpp.Fastq{File: file_handle}
	for {

		_, seq, _, _, err := F.Next()
		seq = bytes.TrimSpace(seq)
		seq = bytes.Replace(seq, []byte("N"), []byte(""), -1)

		for j := 0; j <= len(seq)-kmer; j++ {

			mer := Byte_Int(seq[j : j+kmer])

			total_mer += 1
			_, ok := kmer_number[mer]
			if !ok {
				kmer_spices += 1
				kmer_number[mer] = 1
			} else {
				kmer_number[mer] += 1
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
	average_mer_threshold := total_mer / kmer_spices / 2
	fmt.Println(fmt.Sprintf("Threshold is %d!", average_mer_threshold))
	fmt.Println(len(kmer_number))
	for each_mer, number := range kmer_number {
		if number < average_mer_threshold {
			delete(kmer_number, each_mer)
		}
	}
	fmt.Println(len(kmer_number))
	for e_mer, _ := range kmer_number {
		raw_mer := string(Int_Byte(e_mer, int64(*kmer)))
		rev_mer := lpp.RevComplement([]byte(raw_mer))
		_, ok_rev := kmer_number[Byte_Int(rev_mer)]
		if !ok_rev {
			repre_mer[raw_mer] = ""
		} else {
			all_mer := []string{raw_mer, string(rev_mer)}
			sort.Strings(all_mer)
			repre_mer[all_mer[0]] = ""
		}

		mer_suff := raw_mer[1:*kmer]
		for _, i := range [4]string{"A", "T", "C", "G"} {
			next_mer := mer_suff + i
			_, ok := kmer_number[Byte_Int([]byte(next_mer))]
			if ok {
				_, ok2 := kmer_graph[raw_mer]
				if !ok2 {
					kmer_graph[raw_mer] = make(map[string]int)
				}
				kmer_graph[raw_mer][next_mer] = kmer_number[Byte_Int([]byte(next_mer))]
			}
		}
	}

	for mer1, h1 := range kmer_graph {
		for mer2, number := range h1 {
			Output.WriteString(fmt.Sprintf("%s\t%s\t%d\n", mer1, mer2, number))
		}
	}
}
