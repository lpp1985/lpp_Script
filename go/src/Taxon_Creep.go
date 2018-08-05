package main

import (
	"bufio"
	"flag"
	"fmt"
	. "lpp"
	"os"
)

func Creeper(taxon_number string) {
	all_end[taxon_number] = ""

	_, ok := taxon_data[taxon_number]
	if ok {
		for k1 := range taxon_data[taxon_number] {
			Creeper(k1)
		}
	}
}

func Taxon_Creep(taxon_number string) {
	all_end = make(map[string]string, 1000000)
	Creeper(taxon_number)
}

var all_end map[string]string = make(map[string]string, 1000000)
var taxon_data map[string]map[string]string = make(map[string]map[string]string, 1000000)

func main() {
	defer func() {

		if err := recover(); err != nil {

			switch err {
			case "Input Error":
				fmt.Println(err)
			case "Output Error":
				fmt.Println(err)
			default:
				fmt.Println("aa")
			}
		}
	}()
	file := flag.String("i", "", "input Taxon.dump!")
	taxon := flag.String("n", "", "Taxon Number!")
	output := flag.String("o", "", "Output!")
	flag.Parse()

	Taxon_IO := File_Ddict{File_IO: GetBlockRead(*file, "\n", false, 1000000), Header: false}
	taxon_data = Taxon_IO.Read(3, 1)

	Result_File, err := os.Create(*output)
	if err != nil {
		fmt.Printf("%s is not Exist", *output)
		panic("Output Error")
	}
	defer Result_File.Close()
	BufResult := bufio.NewWriter(Result_File)
	defer BufResult.Flush()
	Taxon_Creep(*taxon)
	for k, _ := range all_end {
		BufResult.WriteString(fmt.Sprintf("%s\t%s\n", k, *taxon))
	}

}
