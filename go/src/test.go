// test.go
package main

import (
	"fmt"
)

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

func Traverse_5(node string, step int, path []string) []string {

	step += 1
	_, ok := kmer_graph[node]
	if !ok || step > 5 {

		return path
	} else {
		path = append(path, node)
		for son, _ := range kmer_graph[node] {
			if Contains(path, son) {
				continue
			}
			fmt.Println(son)

			path = Traverse_5(son, step, path[:step])
			fmt.Println(path)
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
	AddNodes("9", "10")

	fmt.Println(kmer_graph)
	path = Traverse_5("1", 0, path)
	//	fmt.Println(path)
}
