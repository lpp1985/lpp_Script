package main

import (
	"bufio"
	"bytes"
	"flag"
	"lpp"
	"os"
	"strconv"
)

func main() {
	a := flag.String("a", "", "a set")
	b := flag.String("b", "", "b set")
	output := flag.String("o", "", "Output")
	flag.Parse()
	coord_hash := make(map[string]map[int]string)
	RAW := lpp.GetBlockRead(*b, "\n", false, 10000)
	for {
		data, err := RAW.Next()

		data_l := bytes.Split(data, []byte("\t"))
		if len(data_l) > 3 {

			chr := string(data_l[0])
			start, _ := strconv.Atoi(string(data_l[1]))

			end, _ := strconv.Atoi(string(data_l[2]))
			_, ok1 := coord_hash[chr]
			if !ok1 {
				coord_hash[chr] = make(map[int]string)
			}
			for i := start; i <= end; i++ {
				coord_hash[chr][i] = ""

			}
		}
		if err != nil {
			break
		}

	}
	RAW = lpp.GetBlockRead(*a, "\n", false, 10000)
	END, _ := os.Create(*output)
	OUTPUTBUF := bufio.NewWriterSize(END, 10000)
	defer OUTPUTBUF.Flush()
	for {
		data, err := RAW.Next()
		data = bytes.Replace(data, []byte(",<*>"), []byte(""), -1)
		if len(data) > 3 {
			if data[0] == '#' {
				OUTPUTBUF.Write(data)
				continue
			}
			data_l := bytes.Split(data, []byte("\t"))
			chr := string(data_l[0])
			coord, _ := strconv.Atoi(string(data_l[1]))

			_, ok := coord_hash[chr][coord]
			if ok {
				OUTPUTBUF.Write(data)
			}
		}
		if err != nil {
			break
		}
	}
}
