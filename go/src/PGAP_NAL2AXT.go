package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"lpp"
	"os"
	"path/filepath"
	"runtime"
	"strings"
	//	"sync"
)

func Trans(name string, output *string) {
	result_cache := [][]byte{}
	number := strings.Split(filepath.Base(name), ".")[0]
	OUTPUT, err := os.Create(*output + "/" + number + ".axt")
	//	defer OUTPUT.Close()
	if err != nil {
		panic(err)
	}

	BufOUTPUT := bufio.NewWriterSize(OUTPUT, 9999)
	max := 0
	//	defer BufOUTPUT.Flush()
	NAL_HANDLE := lpp.GetBlockRead(name, "\n>", false, 10000)
	BufOUTPUT.WriteString(number + "\n")
	for {

		line, err := NAL_HANDLE.Next()
		if err != nil {
			break
		}
		seq := bytes.SplitN(line, []byte("\n"), 2)[1]
		seq = bytes.TrimSuffix(seq, []byte(">"))
		seq = bytes.Replace(seq, []byte("\n"), []byte(""), -1)
		if bytes.Contains(seq, []byte("N")) {
			continue
		}
		result_cache = append(result_cache, seq)
		for i := 0; i < len(seq); i++ {

			if seq[i] != '-' {
				if i > max {
					max = i
				}
				break
			}
		}

	}
	for _, new_seq := range result_cache {
		if len(new_seq)%2 == 1 {
			new_seq = new_seq[3:]
		}
		BufOUTPUT.Write(new_seq)
		BufOUTPUT.WriteString("\n")
	}
	BufOUTPUT.Flush()
	OUTPUT.Close()
	//	wait.Done()
}
func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())

	//	wg := sync.WaitGroup{}
	defer func() {
		err := recover()

		if err != nil {

			fmt.Println(err)

		}
	}()
	data := flag.String("i", "", "Path")
	output := flag.String("o", "", "Output")
	flag.Parse()

	all_file, err := filepath.Glob(*data + "/*.nal")

	if err != nil {
		panic("nal not exist!!")
	}

	//	wg.Add(len(all_file))

	for _, e_f := range all_file {
		Trans(e_f, output)
	}
	//	wg.Wait()
}
