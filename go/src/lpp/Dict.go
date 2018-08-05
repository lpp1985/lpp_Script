package lpp

import (
	"bytes"
	"sort"
)

type File_Ddict struct {
	File_IO IO
	Header  bool
}

func (file *File_Ddict) Read(key int, value int) map[string]map[string]string {
	key--
	value--
	var result_hash map[string]map[string]string = make(map[string]map[string]string)
	if file.Header == true {
		file.File_IO.Next()
	}

	for {

		line, err := file.File_IO.Next()

		line_l := bytes.Split(bytes.TrimSpace(line), []byte("\t"))

		if len(line_l) >= sort.IntSlice([]int{key, value})[0] {
			key_string := string(line_l[key])
			value_string := string(line_l[value])
			_, ok := result_hash[key_string]

			if !ok {
				result_hash[key_string] = make(map[string]string)

			}
			result_hash[key_string][value_string] = ""
		}
		if err != nil {
			break
		}

	}
	return result_hash
}

type File_dict struct {
	File_IO IO
	Header  bool
}

func (file *File_dict) Read(key int, value int) map[string]string {
	key--
	value--
	var result_hash map[string]string = make(map[string]string)

	if file.Header == true {
		file.File_IO.Next()
	}
	for {

		line, err := file.File_IO.Next()

		line_l := bytes.Split(bytes.TrimSpace(line), []byte("\t"))
		if len(line_l) > sort.IntSlice([]int{key, value})[0] {
			key_string := string(line_l[key])
			value_string := string(line_l[value])
			result_hash[key_string] = value_string

		}
		if err != nil {
			break
		}
	}
	return result_hash
}
