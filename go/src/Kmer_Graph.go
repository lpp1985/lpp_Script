package main

import (
	"bytes"
	"flag"
	"fmt"
	"lpp"

	"os"
	"regexp"
	"strings"
)

var x int

var path []string
var kmer_graph map[string]map[string]string = make(map[string]map[string]string)
var kmer_seq map[string]string = make(map[string]string)

var OUTPUT *os.File

func Contains(s_list []string, node string) bool {
	res := false
	for _, data := range s_list {
		if data == node {
			res = true
		}
	}
	return res
}
func Traverse_5(node string, step int, path []string, threshold int) []string {
	if Contains(path, node) {
		return path
	}
	path = append(path, node)
	_, ok := kmer_graph[node]

	if !ok || step == threshold {
		if len(path) > 0 {
			OUTPUT.WriteString(strings.Join(path, "; ") + "\n")
		}

	} else {
		loc := Find(path, node)
		for son, _ := range kmer_graph[node] {
			if len(kmer_graph[node]) > 1 {
				path = path[:loc]
			}

			path = Traverse_5(son, step, path, threshold)
		}

	}
	return path
}
func Find(s_list []string, node string) int {
	j := 0
	for i, data := range s_list {
		if data == node {
			j = i
			break
		}
	}
	return j + 1
}
func main() {

	input := flag.String("i", "./", "Input File")
	output := flag.String("o", "StringGraph.tsv", "Output File")
	nodes_file := flag.String("l", "KmerList.tsv", "Kmer Id list")
	//	kmer := flag.Int("k", 41, "kmer")
	threshold := flag.Int("s", 10, "Step")
	flag.Parse()
	raw_file := lpp.GetBlockRead(*nodes_file, "\n", false, 10000)
	All_list := new(lpp.File_dict)
	All_list.File_IO = raw_file
	All_list.Header = false
	need_hash := All_list.Read(1, 1)
	for node, _ := range need_hash {
		fmt.Println(node)
	}
	RAW := lpp.Fasta{File: *input}
	OUTPUT, _ = lpp.GetOuput(*output, 1000)

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
	for node, _ := range need_hash {
		OUTPUT.WriteString(node + "\n")
		Traverse_5(node, 0, path, *threshold)

	}
}
