package main

import (
	"bytes"
	"lpp"
	"os"
	"regexp"
)

func main() {
	FASTAIO := lpp.GetBlockRead(os.Args[1], "\n>", false, 10000)
	OUTPUTIO, _ := lpp.GetOuput(os.Args[2], 100000)
	reg := regexp.MustCompile(`(H\wH\wDH)`)
	for {
		block_data, err := FASTAIO.Next()
		block_data = bytes.TrimSuffix(block_data, []byte(">"))

		if err != nil {
			break
		}
		match := bytes.SplitN(block_data, []byte("\n"), 2)
		seq := match[1]
		seq2 := bytes.Replace(seq, []byte("\n"), []byte(""), -1)
		seq2 = bytes.ToUpper(seq2)
		if bytes.Contains(seq2, []byte("X")) {
			continue
		}
		result := reg.Find(seq2)
		if len(result) > 0 {
			if string(block_data[0]) != ">" {
				block_data = append([]byte(">"), block_data...)
			}

			OUTPUTIO.Write(block_data)

		}
	}

}
