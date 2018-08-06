// test.go
package main

import (
	"bytes"
	"fmt"
	"lpp"
	"os"
	"strconv"
)

func main() {
	RAWHANDLE, _ := os.Open(os.Args[1])
	RAWIO := lpp.GetBlockRead(RAWHANDLE, "\n", false, 10000)
	max := 0
	all_need := make(map[string]string, 10000)
	all_cov := 0
	END, _ := os.Create(os.Args[2])
	raw_name := ""
	for {
		line, err := RAWIO.Next()

		line_l := bytes.Fields(line)
		c_name := string(line_l[0])
		_, ok := all_need[c_name]
		//fmt.Println(ok)

		if ok == false || err != nil {

			if max != 0 && raw_name != "" {
				fmt.Println(ok, raw_name, err)
				t_cov := all_cov / max
				END.WriteString(fmt.Sprintf("%s\t%d\n", raw_name, t_cov))

				max = 0
				all_cov = 0
			}
			raw_name = c_name
		}
		if err != nil {
			break
		}
		all_need[c_name] = ""
		coord, _ := strconv.Atoi(string(line_l[1]))
		cov, err := strconv.Atoi(string(line_l[2]))
		fmt.Println(cov)
		max = coord
		all_cov += cov

	}
}
