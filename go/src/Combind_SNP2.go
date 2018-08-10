package main

import (
	"bufio"
	"bytes"
	"flag"
	//	"fmt"
	//	"io/ioutil"
	"lpp"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
)

func main() {

	all_data, _ := filepath.Glob("*.v2f")
	output := flag.String("o", "", "Output")
	flag.Parse()
	coord_hash := make(map[string]map[int]map[string]string)
	all_sample := make([]string, 0, 100)
	all_sample = append(all_sample, "Ref")
	var all_fasta []*bytes.Buffer = make([]*bytes.Buffer, 1, 1000)
	//	data, _ := os.Create("Ref.fa")
	//	CACHEBUF := bufio.(data, 1000000)
	//	CACHEBUF.WriteString(">" + "Ref" + "\n")
	//	defer CACHEBUF.Flush()

	all_fasta[0] = new(bytes.Buffer)
	all_fasta[0].WriteString(">" + "Ref" + "\n")
	for _, eachfile := range all_data {
		RAW := lpp.GetBlockRead(eachfile, "\n", false, 10000)
		name := strings.Split(eachfile, ".")[0]
		all_sample = append(all_sample, name)

		cache := new(bytes.Buffer)
		cache.WriteString(">" + "Ref" + "\n")
		//		CACHEBUF := bufio.NewWriterSize(data, 1000000)
		//		CACHEBUF.WriteString(">" + name + "\n")
		//		defer CACHEBUF.Flush()
		all_fasta = append(all_fasta, cache)

		for {
			data, err := RAW.Next()

			data_l := bytes.Split(data, []byte("\t"))
			if len(data_l) > 3 {
				if data[0] == '#' {
					continue
				}
				chr := string(data_l[0])
				coord, _ := strconv.Atoi(string(data_l[1]))

				_, ok1 := coord_hash[chr]
				if !ok1 {
					coord_hash[chr] = make(map[int]map[string]string)
				}
				_, ok2 := coord_hash[chr][coord]
				if !ok2 {
					coord_hash[chr][coord] = make(map[string]string)
				}
				ref_site := string(data_l[3])
				alt_site := string(bytes.Split(data_l[4], []byte(","))[0])
				coord_hash[chr][coord]["Ref"] = ref_site
				coord_hash[chr][coord][name] = alt_site
			}
			if err != nil {
				break
			}

		}
	}

	END, _ := os.Create(*output)
	OUTPUTBUF := bufio.NewWriterSize(END, 10000)
	defer OUTPUTBUF.Flush()
	chr_all := make([]string, 0, 1000000)
	for e_chr, _ := range coord_hash {
		chr_all = append(chr_all, e_chr)
	}
	sort.Strings(chr_all)

	for _, each_chr := range chr_all {
		site_all := make([]int, 0, 1000000000)
		site_hash, _ := coord_hash[each_chr]

		for e_coord, _ := range site_hash {
			site_all = append(site_all, e_coord)
		}

		sort.Ints(site_all)
		for _, site := range site_all {
			for i, e_sample := range all_sample {

				char, data_ok := site_hash[site][e_sample]

				if !data_ok {
					char = site_hash[site]["Ref"]
				}

				all_fasta[i].WriteString(char)
			}
		}

	}
	for _, handle := range all_fasta {
		OUTPUTBUF.WriteString(handle.String() + "\n")

	}

}
