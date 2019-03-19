package main

import (
	//"fmt"
	"bytes"
	"flag"
	"lpp"
	"regexp"
	"strconv"

	mgo "gopkg.in/mgo.v2"
	"gopkg.in/mgo.v2/bson"
)

var kmer_graph map[string]map[string]string = make(map[string]map[string]string)

type Kmer struct {
	Name string
	Son  map[string]string
}

type Sequence struct {
	Name   string
	ForSeq string
	Rever  string
}

var session, err = mgo.Dial("192.168.31.82:27017")

func Insert(table *mgo.Collection, query string, subject string) {
	result := Kmer{}
	errfind := table.Find(bson.M{"name": query}).One(&result)
	if errfind != nil {
		table.Insert(&Kmer{Name: query, Son: map[string]string{subject: ""}})

	} else {
		table.Find(bson.M{"name": query}).One(&result)
		result.Son[subject] = ""
		err := table.Update(bson.M{"name": query}, result)
		if err != nil {
			panic(err)
		}
	}

}
func main() {
	has_son := make(map[string]string)
	fasta := flag.String("i", "./", "Input String Fasta")
	name := flag.String("n", "Kmer", "DB Name")
	flag.Parse()
	session.SetMode(mgo.Monotonic, true)
	defer session.Close()
	var db_kmer = session.DB(*name).C("Kmer")
	var db_seq = session.DB(*name).C("Sequence")
	db_kmer.RemoveAll(nil)
	db_seq.RemoveAll(nil)
	RAW := lpp.Fasta{File: *fasta}
	reg := regexp.MustCompile(`L\:(\S)\:(\d+)\:(\S)`)
	reg_cov := regexp.MustCompile(`KC\:i\:(\d+)`)
	all_seq := make([]interface{}, 0, 1000000)
	x := 0
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
		all_seq = append(all_seq, &Sequence{Name: name, ForSeq: string(seq), Rever: string(lpp.RevComplement(seq))})
		x += 1
		//fmt.Println( "X is ", x    )
		if x == 1000000 {
			//fmt.Println( len(all_seq) )
			//fmt.Println(x)
			x = 0
			db_seq.Insert(all_seq...)
			all_seq = make([]interface{}, 0, 1000000)
			//fmt.Println( len(all_seq) )
		}
		if err != nil {
			break
		}
	}
	if len(all_seq) > 0 {
		db_seq.Insert(all_seq...)
	}
	RAW = lpp.Fasta{File: *fasta}
	for {
		title, _, err := RAW.Next()
		if len(title) == 0 {
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
	for query, hash := range kmer_graph {
		db_kmer.Insert(&Kmer{Name: query, Son: hash})
	}

}
