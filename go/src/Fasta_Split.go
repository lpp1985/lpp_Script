// Fasta_Split.go
package main

import (
	"lpp"
	"os"
)

func main() {
	data := lpp.Fasta{File: os.Args[0]}
	END, _ := lpp.GetOuput(os.Args[1], 10000000)
	i := 0
	for {
		i += 1
		t, s, err := data.Next()
		END.Write(t)
		END.Write(s)
		if i > 250000 {
			break
		}
		if err != nil {
			break
		}
	}
}
