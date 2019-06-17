<<<<<<< HEAD
package main

import (
	"bytes"
	"flag"

	//"fmt"
	"lpp"
	"os"
	"regexp"
	"strconv"
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
func Traverse_5(order int, node string, step int, path []string, threshold int, OUT *os.File) []string {
	if Contains(path, node) {
		return path
	}
	path = append(path, node)
	//	_, ok := kmer_graph[node]

	if step == threshold {
		status := 0
		start := path[0][:len(path)-1]
		for _, name := range path[1:] {
			q_name := name[:len(name)-1]
			if q_name == start {
				status = 1
			}

		}
		if len(path) > 0 && status == 0 {
			OUTPUT.WriteString(">" + strings.Join(path, "; ") + "\n")
			end_seq := ""
			for i, name := range path {
				new_seq, _ := kmer_seq[name]
				if i == 0 {

					end_seq = new_seq
				} else {
					end_seq += new_seq[*kmer-1:]
				}
			}
			if order%2 == 1 {
				end_seq = string(lpp.RevComplement([]byte(end_seq)))
			}
			OUTPUT.WriteString(end_seq + "\n")
		}

	} else {
		loc := Find(path, node)
		if len(kmer_graph[node]) > 1 {
			step += 1
		}
		for son, _ := range kmer_graph[node] {
			if len(kmer_graph[node]) > 1 {
				path = path[:loc]
			}

			path = Traverse_5(order, son, step, path, threshold, OUT)
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

var kmer *int

func main() {
	has_son := make(map[string]string)
	input := flag.String("i", "./", "Input File")
	output_forward := flag.String("f", "StringGraph_for.tsv", "Output File")
	output_reverse := flag.String("r", "StringGraph_rev.tsv", "Output File")
	nodesforward_file := flag.String("l", "KmerList_for.tsv", "Kmer Id list")
	nodesreverse_file := flag.String("j", "KmerList_rev.tsv", "Kmer Id list")
	kmer = flag.Int("k", 41, "kmer")
	threshold := flag.Int("s", 10, "Step")
	flag.Parse()
	raw_for_file := lpp.GetBlockRead(*nodesforward_file, "\n", false, 10000)
	raw_rev_file := lpp.GetBlockRead(*nodesreverse_file, "\n", false, 10000)
	var All_for_list []string
	var All_rev_list []string
	for {
		line, err := raw_for_file.Next()
		line = bytes.TrimSpace(line)
		All_for_list = append(All_for_list, string(line))
		if err != nil {
			break
		}

	}
	for {
		line, err := raw_rev_file.Next()
		line = bytes.TrimSpace(line)
		All_rev_list = append(All_rev_list, string(line))
		if err != nil {
			break
		}

	}
	RAW := lpp.Fasta{File: *input}
	OUTPUT_FOR, _ := lpp.GetOuput(*output_forward, 1000)
	OUTPUT_REV, _ := lpp.GetOuput(*output_reverse, 1000)
	reg := regexp.MustCompile(`L\:(\S)\:(\d+)\:(\S)`)
	reg_cov := regexp.MustCompile(`KC\:i\:(\d+)`)

	for {
		title, seq, err := RAW.Next()
		if len(title) == 0 {
			break
		}
		name := string(bytes.Fields(title)[0][1:])
		all_situation := reg.FindAllStringSubmatch(string(title), -1)
		cov, _ := strconv.Atoi(reg_cov.FindStringSubmatch(string(title))[1])
		if cov < 30 {
			continue
		}
		if len(all_situation) == 0 {
			continue
		}
		if len(all_situation) < 2 {
			continue
		}
		has_son[name] = ""
		seq = bytes.TrimSpace(seq)
		kmer_seq[name+"+"] = string(seq)
		kmer_seq[name+"-"] = string(lpp.RevComplement(seq))
		if err != nil {
			break
		}
	}
	RAW = lpp.Fasta{File: *input}
	for {
		title, _, err := RAW.Next()
		name := string(bytes.Fields(title)[0][1:])
		cov, _ := strconv.Atoi(reg_cov.FindStringSubmatch(string(title))[1])
		if cov < 30 {
			continue
		}
		_, ok := has_son[name]
		if !ok {
			continue
		}
		all_situation := reg.FindAllStringSubmatch(string(title), -1)
		for _, data := range all_situation {
			dir := data[1]
			sub := data[2]
			_, ok := has_son[sub]
			if !ok {
				continue
			}
			sub_dir := data[3]
			q := name + dir
			s := sub + sub_dir
			_, ok1 := kmer_graph[q]

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

		}
		if err != nil {
			break
		}
	}
	for _, node := range All_for_list {
		// OUTPUT.WriteString(node + "\n")
		Traverse_5(0, node, 0, path, *threshold, OUTPUT_FOR)

	}
	for _, node := range All_rev_list {
		// OUTPUT.WriteString(node + "\n")
		Traverse_5(0, node, 0, path, *threshold, OUTPUT_REV)

	}
}
=======
package main

import (
	"bytes"
	"flag"

	"fmt"
	"lpp"
	"os"
	"regexp"
	"strconv"
	"strings"
)

var x int

var path []string
var kmer_graph map[string]map[string]string = make(map[string]map[string]string)
var kmer_seq map[string]string = make(map[string]string)


func Contains(s_list []string, node string) bool {
	res := false
	for _, data := range s_list {
		if data == node {
			res = true
		}
	}
	return res
}
func Traverse_5(order int, node string, step int, path []string, threshold int, OUT *os.File) []string {
	if Contains(path, node) {
		return path
	}
	path = append(path, node)
	//	_, ok := kmer_graph[node]

	if step == threshold {
		status := 0
		
		start := path[0][:len(path[0])-1]
		for _, name := range path[1:] {
			q_name := name[:len(name)-1]
			if q_name == start {
				status = 1
			}

		}
		fmt.Println(path,status)
		if len(path) > 0 && status == 0 {
			
			OUT.WriteString(">" + strings.Join(path, "; ") + "\n")
			end_seq := ""
			for i, name := range path {
				new_seq, _ := kmer_seq[name]
				if i == 0 {

					end_seq = new_seq
				} else {
					end_seq += new_seq[*kmer-1:]
				}
			}
			if order%2 == 1 {
				end_seq = string(lpp.RevComplement([]byte(end_seq)))
			}
			OUT.WriteString(end_seq + "\n")
		}

	} else {
		loc := Find(path, node)
		if len(kmer_graph[node]) > 1 {
			step += 1
		}
		for son, _ := range kmer_graph[node] {
			if len(kmer_graph[node]) > 1 {
				path = path[:loc]
			}

			path = Traverse_5(order, son, step, path, threshold, OUT)
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

var kmer *int

func main() {
	has_son := make(map[string]string)
	input := flag.String("i", "./", "Input File")
	output_forward := flag.String("f", "StringGraph_for.tsv", "Output File")
	output_reverse := flag.String("r", "StringGraph_rev.tsv", "Output File")
	nodesforward_file := flag.String("l", "KmerList_for.tsv", "Kmer Id list")
	nodesreverse_file := flag.String("j", "KmerList_rev.tsv", "Kmer Id list")
	kmer = flag.Int("k", 41, "kmer")
	threshold := flag.Int("s", 10, "Step")
	flag.Parse()
	raw_for_file := lpp.GetBlockRead(*nodesforward_file, "\n", false, 10000)
	raw_rev_file := lpp.GetBlockRead(*nodesreverse_file, "\n", false, 10000)
	var All_for_list []string
	var All_rev_list []string
	for {
		line, err := raw_for_file.Next()
		line = bytes.TrimSpace(line)
		All_for_list = append(All_for_list, string(line))
		if err != nil {
			break
		}

	}
	for {
		line, err := raw_rev_file.Next()
		line = bytes.TrimSpace(line)
		All_rev_list = append(All_rev_list, string(line))
		if err != nil {
			break
		}

	}
	RAW := lpp.Fasta{File: *input}
	OUTPUT_FOR, _ := lpp.GetOuput(*output_forward, 1000)
	OUTPUT_REV, _ := lpp.GetOuput(*output_reverse, 1000)
	reg := regexp.MustCompile(`L\:(\S)\:(\d+)\:(\S)`)
	reg_cov := regexp.MustCompile(`KC\:i\:(\d+)`)

	for {
		title, seq, err := RAW.Next()
		if len(title) == 0 {
			break
		}
		name := string(bytes.Fields(title)[0][1:])
		all_situation := reg.FindAllStringSubmatch(string(title), -1)
		cov, _ := strconv.Atoi(reg_cov.FindStringSubmatch(string(title))[1])
		if cov < 30 {
			continue
		}
		if len(all_situation) == 0 {
			continue
		}
		if len(all_situation) < 2 {
			continue
		}
		has_son[name] = ""
		seq = bytes.TrimSpace(seq)
		kmer_seq[name+"+"] = string(seq)
		kmer_seq[name+"-"] = string(lpp.RevComplement(seq))
		if err != nil {
			break
		}
	}
	RAW = lpp.Fasta{File: *input}
	for {
		title, _, err := RAW.Next()
		if len(title)<1{
			break
		}
		name := string(bytes.Fields(title)[0][1:])
		cov, _ := strconv.Atoi(reg_cov.FindStringSubmatch(string(title))[1])
		if cov < 30 {
			continue
		}
		_, ok := has_son[name]
		if !ok {
			continue
		}
		all_situation := reg.FindAllStringSubmatch(string(title), -1)
		for _, data := range all_situation {
			dir := data[1]
			sub := data[2]
			_, ok := has_son[sub]
			if !ok {
				continue
			}
			sub_dir := data[3]
			q := name + dir
			s := sub + sub_dir
			_, ok1 := kmer_graph[q]

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

		}
		if err != nil {
			break
		}
	}
	for _, node := range All_for_list {
		// OUTPUT.WriteString(node + "\n")
		Traverse_5(0, node, 0, path, *threshold, OUTPUT_FOR)

	}
	for _, node := range All_rev_list {
		// OUTPUT.WriteString(node + "\n")
		Traverse_5(0, node, 0, path, *threshold, OUTPUT_REV)

	}
}
>>>>>>> bd03fdf8f52fee002da8159f1b0eb48ece0d0468
