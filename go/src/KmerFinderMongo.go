// KmerFinder
package main

import (
	"bytes"
	"flag"
	"fmt"
	"lpp"

	"runtime"
	"sync"

	"gopkg.in/mgo.v2"
	"gopkg.in/mgo.v2/bson"
)

var session, _ = mgo.Dial("192.168.31.82:27017")
var db = session.DB("Kmer").C("Kmer")

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

//var kmer_number map[string]int = make(map[string]int)
var kmer_graph map[string]map[string]int = make(map[string]map[string]int)

func MultiProcess(wg *sync.WaitGroup, seq_list []interface{}) {

	db.Insert(seq_list...)
	wg.Done()

}
func Kmer_Count(file_handle string) {
	wg := sync.WaitGroup{}
	wg.Add(64)
	F := lpp.Fasta{File: file_handle}
	multiprocess_hash := make([][]interface{}, 0, 64)
	m := 10000
	kmer_cache := make([]interface{}, 0, m)
	i := 1
	j := 1

	for {

		_, seq, err := F.Next()
		seq = bytes.TrimSpace(seq)
		is_n := bytes.Contains(seq, []byte("N"))
		if is_n {
			continue
		}

		kmer_cache = append(kmer_cache, &KMER{Kmer: string(seq)})

		if i == m {
			multiprocess_hash = append(multiprocess_hash, kmer_cache)

			kmer_cache = make([]interface{}, 0, m)
			i = 0
			if j == 64*m {
				fmt.Println(len(multiprocess_hash))
				for _, data := range multiprocess_hash {
					go MultiProcess(&wg, data)

				}

				wg.Wait()
				multiprocess_hash = make([][]interface{}, 0, 64)
				wg.Add(64)
				j = 0
			}

		}
		i += 1
		j += 1

		if err != nil {

			break
		}
	}
	for _, data := range multiprocess_hash {
		db.Insert(data...)
	}
}

type KMER struct {
	Kmer string
}

func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())
	read1 := flag.String("f", "", "read1")
	kmer := flag.Int("k", 41, "kmer")
	output := flag.String("o", "", "Output")

	flag.Parse()
	index := mgo.Index{Key: []string{"$text:kmer"}}

	err := db.EnsureIndex(index)
	if err != nil {
		panic(err)
	}
	session.SetMode(mgo.Monotonic, true)
	defer session.Close()
	Output, _ := lpp.GetOuput(*output, 10000)
	fmt.Println("read1 Start!!")
	Kmer_Count(*read1)
	kmer_number := make([]KMER, 0, 10000000000)
	db.Find(nil).All(&kmer_number)
	fmt.Println(len(kmer_number))
	for _, e_mer := range kmer_number {
		raw_mer := e_mer.Kmer
		mer_suff := raw_mer[1:*kmer]
		for _, i := range [4]string{"A", "T", "C", "G"} {
			next_mer := mer_suff + i
			result := KMER{}
			err := db.Find(bson.M{"kmer": next_mer}).One(&result)

			if err == nil {
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

//	for raw_mer, _ := range kmer_number {

//		mer_suff := raw_mer[1:*kmer]
//		for _, i := range [4]string{"A", "T", "C", "G"} {
//			next_mer := mer_suff + i
//			_, ok := kmer_number[next_mer]
//			if ok {
//				_, ok2 := kmer_graph[raw_mer]
//				if !ok2 {
//					kmer_graph[raw_mer] = make(map[string]int)
//				}
//				kmer_graph[raw_mer][next_mer] = 0
//			}
//		}
//	}

//	for mer1, h1 := range kmer_graph {
//		for mer2, number := range h1 {
//			Output.WriteString(fmt.Sprintf("%s\t%s\t%d\n", mer1, mer2, number))
//		}
//	}
