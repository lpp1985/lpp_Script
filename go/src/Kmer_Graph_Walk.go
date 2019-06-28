package main

import (
	"bytes"
	"flag"

	//"fmt"

	"lpp"
	"os"
	"strings"

	mgo "gopkg.in/mgo.v2"
	"gopkg.in/mgo.v2/bson"
)

type Kmer struct {
	Name string
	Son  map[string]string
}

var session, _ = mgo.Dial("192.168.31.82:27017")

var path []string

var OUTPUT *os.File

type Sequence struct {
	Name string
	Seq  string
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

		if len(path) > 0 && status == 0 {

			OUT.WriteString(">" + strings.Join(path, "; ") + "\n")
			end_seq := ""
			for i, name := range path {
				seq_searchresult := Sequence{}
				db_seq.Find(bson.M{"name": name}).One(&seq_searchresult)

				new_seq := seq_searchresult.Seq
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

		search_result := Kmer{}
		db_kmer.Find(bson.M{"name": node}).One(&search_result)
		if len(search_result.Son) > 1 {
			step += 1
		}
		for son, _ := range search_result.Son {
			if len(search_result.Son) > 1 {
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

// func Traverse_5(order int, node string, step int, path []string, threshold int, OUT *os.File) []string {
// 	if Contains(path, node) {
// 		return path
// 	}
// 	path = append(path, node)
// 	//	_, ok := kmer_graph[node]

// 	if step == threshold {
// 		if len(path) > 0 {
// 			OUTPUT.WriteString(">" + strings.Join(path, "; ") + "\n")
// 			end_seq := ""
// 			for i, name := range path {
// 				seq_searchresult := Sequence{}
// 				db_seq.Find(bson.M{"name": name}).One(&seq_searchresult)

// 				new_seq := seq_searchresult.Seq
// 				if i == 0 {

// 					end_seq = new_seq
// 				} else {
// 					end_seq += new_seq[*kmer-1:]
// 				}
// 			}
// 			if order%2 == 1 {
// 				end_seq = string(lpp.RevComplement([]byte(end_seq)))
// 			}
// 			OUTPUT.WriteString(end_seq + "\n")
// 		}

// 	} else {
// 		loc := Find(path, node)

// 		search_result := Kmer{}
// 		db_kmer.Find(bson.M{"name": node}).One(&search_result)
// 		if len(search_result.Son) > 1 {
// 			step += 1
// 		}
// 		for son, _ := range search_result.Son {
// 			if len(search_result.Son) > 1 {
// 				path = path[:loc]
// 			}

// 			path = Traverse_5(order, son, step, path, threshold)
// 		}

// 	}
// 	return path
// }

var kmer *int

var name *string
var db_kmer *mgo.Collection
var db_seq *mgo.Collection

func main() {

	name = flag.String("i", "Sample", "Sample Name")

	output_forward := flag.String("f", "StringGraph_for.tsv", "Output File")
	output_reverse := flag.String("r", "StringGraph_rev.tsv", "Output File")
	nodesforward_file := flag.String("l", "KmerList_for.tsv", "Kmer Id list")
	nodesreverse_file := flag.String("j", "KmerList_rev.tsv", "Kmer Id list")
	kmer = flag.Int("k", 41, "kmer")
	threshold := flag.Int("s", 10, "Step")
	flag.Parse()
	db_kmer = session.DB(*name).C("Kmer")
	db_seq = session.DB(*name).C("Sequence")
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

	OUTPUT_FOR, _ := lpp.GetOuput(*output_forward, 1000)
	OUTPUT_REV, _ := lpp.GetOuput(*output_reverse, 1000)

	for _, node := range All_for_list {
		// OUTPUT.WriteString(node + "\n")
		Traverse_5(0, node, 0, path, *threshold, OUTPUT_FOR)

	}
	for _, node := range All_rev_list {
		// OUTPUT.WriteString(node + "\n")
		Traverse_5(0, node, 0, path, *threshold, OUTPUT_REV)
	}
}
