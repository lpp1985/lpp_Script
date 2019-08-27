package main

import (
	"bytes"
	"flag"

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

type Node struct {
	Name string
}
type Start struct {
	Name string
	Attr string
}

var AllrelatedNodes map[string]string
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
			path1 := []string{}
			for i := 0; i < len(path); i++ {
				path1 = append(path1, Rename(path[i]))
			}
			OUT.WriteString(">" + strings.Join(path, "; ") + "\n")
			end_seq := []string{}

			for i, name := range path {
				db_related.Insert(Start{Name: name, Attr: path[0]})
				seq_searchresult := Sequence{}
				db_seq.Find(bson.M{"name": name}).One(&seq_searchresult)

				new_seq := seq_searchresult.Seq
				if i == 0 {

					end_seq = append(end_seq, new_seq)
				} else {
					end_seq = append(end_seq, new_seq[*kmer-1:])
				}
			}
			res_order := []string{}
			if order%2 == 1 {

				for i := len(end_seq) - 1; i >= 0; i-- {
					res_order = append(res_order, string(lpp.RevComplement([]byte(end_seq[i]))))
				}
				k := 0
				for i := len(end_seq) - 1; i >= 0; i-- {
					res_order[k] = res_order[k] + "(" + path[i] + ")"
					k += 1
				}

			} else {
				for i := 0; i < len(end_seq); i++ {
					res_order = append(res_order, end_seq[i])
				}
				k := 0
				for i := 0; i < len(end_seq); i++ {
					res_order[k] = res_order[k] + "(" + path[i] + ")"
					k += 1
				}
			}

			OUT.WriteString(strings.Join(res_order, "---------") + "\n")
		}

	} else {
		loc := Find(path, node)

		search_result := Kmer{}

		db_kmer.Find(bson.M{"name": node}).One(&search_result)
		has := search_result.Son
		for each_data, _ := range has {
			s_res := Node{}
			err := db_delete.Find(bson.M{"name": each_data}).One(&s_res)
			if err == nil {
				delete(has, each_data)

			}
		}
		if len(has) > 1 {
			step += 1
		}
		for son, _ := range has {
			if len(has) > 1 {
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

var name *string
var db_kmer *mgo.Collection
var db_seq *mgo.Collection
var db_start *mgo.Collection
var db_delete *mgo.Collection
var db_related *mgo.Collection

func Rename(x string) string {
	data := []byte(x)
	length := len(data) - 1
	if data[length] == '+' {
		data[length] = '-'
	} else {
		data[length] = '-'
	}
	return string(data)
}
func main() {
	index := mgo.Index{
		Key:      []string{"name", "attr"},
		Unique:   true,
		DropDups: true,
	}
	name = flag.String("i", "Sample", "Sample Name")

	output_forward := flag.String("f", "StringGraph_for.tsv", "Output File")
	output_reverse := flag.String("r", "StringGraph_rev.tsv", "Output File")
	nodesforward_file := flag.String("l", "", "Kmer Id list")
	nodesreverse_file := flag.String("j", "", "Kmer Id list")
	flush := flag.Bool("b", false, "Is flush the database")
	kmer = flag.Int("k", 41, "kmer")
	threshold := flag.Int("s", 10, "Step")
	flag.Parse()
	db_kmer = session.DB(*name).C("Kmer")
	db_related = session.DB(*name).C("Related")
	db_start = session.DB(*name).C("Start")
	db_delete = session.DB(*name).C("Delete")
	db_delete.EnsureIndex(index)
	db_start.EnsureIndex(index)
	db_related.EnsureIndex(index)
	db_seq = session.DB(*name).C("Sequence")
	raw_for_file := lpp.IO{}
	raw_rev_file := lpp.IO{}
	var All_for_list []string
	var All_rev_list []string
	if *nodesforward_file != "" && *nodesreverse_file != "" {
		if *flush {
			db_start.RemoveAll(nil)
		}
		raw_for_file = lpp.GetBlockRead(*nodesforward_file, "\n", false, 10000)
		raw_rev_file = lpp.GetBlockRead(*nodesreverse_file, "\n", false, 10000)
		for {
			line, err := raw_for_file.Next()
			line = bytes.TrimSpace(line)
			if string(line) != "" {
				All_for_list = append(All_for_list, string(line))
			}

			if err != nil {
				break
			}

		}
		for _, node := range All_for_list {
			if *flush {

				db_start.Insert(Start{Name: node, Attr: "5'"})
			}
		}
		for {
			line, err := raw_rev_file.Next()
			line = bytes.TrimSpace(line)
			if string(line) != "" {
				All_rev_list = append(All_rev_list, string(line))
			}

			if err != nil {
				break
			}

		}
		for _, node := range All_rev_list {
			if *flush {
				db_start.Insert(Start{Name: node, Attr: "3'"})
			}

		}
	} else {
		start_node := Start{}
		iter := db_start.Find(nil).Iter()

		for iter.Next(&start_node) {

			if start_node.Attr == "5'" {
				All_for_list = append(All_for_list, string(start_node.Name))
			} else {
				All_rev_list = append(All_rev_list, string(start_node.Name))
			}
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
		Traverse_5(1, node, 0, path, *threshold, OUTPUT_REV)
	}
}
