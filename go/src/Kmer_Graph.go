package main

import (
	"bytes"
	//	"encoding/binary"
	//	"fmt"
	"lpp"
	//	"math"

	"os"
	"regexp"
	//	"strconv"
)

var kmer_graph map[string]map[string]string = make(map[string]map[string]string)
var kmer_seq map[string]string = make(map[string]string)

func main() {

	RAW := lpp.Fasta{File: os.Args[1]}
	reg := regexp.MustCompile(`L\:(\S)\:(\d+)\:(\S)`)

	for {
		title, seq, err := RAW.Next()
		name := string(bytes.Fields(title)[0][1:])
		all_situation := reg.FindAllStringSubmatch(string(title), -1)
		if len(all_situation) == 0 {
			continue
		}
		for _, data := range all_situation {
			dir := data[1]
			sub := data[2]
			sub_dir := data[3]
			q := name + dir
			s := sub + sub_dir
			_, ok1 := kmer_graph[q]
			//			fmt.Println(q, s)

			if !ok1 {
				kmer_graph[q] = make(map[string]string)

			}
			kmer_graph[q][s] = ""
			if dir == "+" {
				dir = "-"
			} else {
				dir = "+"
			}
			if sub_dir == "+" {
				sub_dir = "-"
			} else {
				sub_dir = "+"
			}
			q = sub + sub_dir
			s = name + dir
			_, ok2 := kmer_graph[q]
			if !ok2 {
				kmer_graph[q] = make(map[string]string)

			}
			kmer_graph[q][s] = ""
			//			fmt.Println(q, s)
		}

		seq = bytes.TrimSpace(seq)
		kmer_seq[name+"+"] = string(seq)
		kmer_seq[name+"-"] = string(lpp.RevComplement(seq))

		if err != nil {
			break
		}
	}
}
func Contains(s_list []string, node string) bool {
	res := false
	for _, data := range s_list {
		if data == node {
			res = true
		}
	}
	return res
}
func Traverse_5(node string, step int, path []string) []string {

	path = append(path, node)
	step += 1
	_, ok := kmer_graph[node]
	if !ok || step == 5 || Contains(path, node) {
		return path
	} else {
		for _, son := range kmer_graph[node] {
			Traverse_5(son, step, path)
		}

	}
	return path
}
