// test.go
package main

import (
	"fmt"
)

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
func Contains(s_list []string, node string) bool {
	res := false
	for _, data := range s_list {
		if data == node {
			res = true
		}
	}
	return res
}

var kmer_graph map[string]map[string]string = make(map[string]map[string]string)

var x int

func Traverse_5(node string, step int, path []string) []string {
	if Contains(path, node) {
		return path
	}

	path = append(path, node)
	//	_, ok := kmer_graph[node]
	//	fmt.Println(node, step)
	if step == 2 {
		if len(path) > 0 {
			fmt.Println(path, step)
			return path
		}

	} else {
		if len(kmer_graph[node]) > 1 {
			step += 1
		}
		for son, _ := range kmer_graph[node] {
			loc := Find(path, node)
			if len(kmer_graph[node]) > 1 {

				//				fmt.Println(path, loc, node, son)
				path = path[:loc]
			}

			path = Traverse_5(son, step, path)
		}

	}
	return path
}

var path []string

func AddNodes(n1 string, n2 string) {
	_, ok := kmer_graph[n1]
	if !ok {
		kmer_graph[n1] = make(map[string]string)
	}
	kmer_graph[n1][n2] = ""
}

func main() {
	aa := []string{"1", "2", "3", "4", "5", "6", "7"}
	//	aa := []string{"1", "2", "3", "4"}
	for i, _ := range aa {
		if i == len(aa)-1 {
			break
		}
		//		fmt.Println(aa[i], aa[i+1])
		AddNodes(aa[i], aa[i+1])

	}
	AddNodes("3", "2")
	AddNodes("3", "8")
	AddNodes("8", "9")
	AddNodes("9", "11")
	AddNodes("9", "10")
	AddNodes("3", "12")
	AddNodes("12", "13")
	AddNodes("12", "14")
	AddNodes("13", "15")
	AddNodes("14", "15")
	AddNodes("15", "16")
	AddNodes("16", "17")
	AddNodes("15", "17")

	//	fmt.Println(kmer_graph)
	path = Traverse_5("1", 0, path)

}
