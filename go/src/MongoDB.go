package main

import (
	//	"log"
	"flag"
	"fmt"
	"os/exec"
	//	"bufio"
	"bytes"
	//	"crypto/md5"
	//	"io"
	"os"

	"path/filepath"
	"runtime"
	"sync"

	"gopkg.in/mgo.v2"
	"gopkg.in/mgo.v2/bson"
)

type File struct {
	Name string
	Md5  string
}
type File_d struct {
	Name string
}

//func md5sum3(file string) string {
//	f, err := os.Open(file)
//	if err != nil {
//		return ""
//	}
//	defer f.Close()
//	r := bufio.NewReaderSize(f, 10^6)

//	h := md5.New()

//	_, err = io.Copy(h, r)
//	if err != nil {
//		return ""
//	}
//	return fmt.Sprintf("%x", h.Sum(nil))

//}

func getFilelist(path string) []string {

	file_list := []string{}
	err := filepath.Walk(path, func(path string, f os.FileInfo, err error) error {
		if f == nil {
			return err
		}
		if f.IsDir() {
			return nil
		}
		f_path, _ := filepath.Abs(path)
		file_list = append(file_list, f_path)

		return nil
	})
	if err != nil {
		fmt.Printf("filepath.Walk() returned %v\n", err)
	}
	return file_list

}

func Md5_Check(wg *sync.WaitGroup, file string) {
	cmd := exec.Command("md5sum", "-b", fmt.Sprintf("%s", file))
	out, err := cmd.Output()
	if err == nil {
		md5 := string(bytes.Fields(out)[0])
		result := File{}
		errfind := db.Find(bson.M{"md5": md5}).One(&result)
		if errfind != nil {
			db.Insert(&File{Name: file, Md5: md5})

		} else {
			errdel_find := db_delete.Find(bson.M{"name": file}).One(&result)
			if errdel_find != nil {
				db_delete.Insert(&File_d{Name: file})
			}

		}

	} else {
		fmt.Println(err)
	}
	wg.Done()
}

var session, err = mgo.Dial("192.168.0.106:27017")

var db = session.DB("MD5").C("AV")

var db_delete = session.DB("MD5").C("Delete")

func main() {
	runtime.GOMAXPROCS(runtime.NumCPU())
	wg := sync.WaitGroup{}
	start_path := flag.String("i", "./", "Input Path")
	flag.Parse()
	defer session.Close()
	if err != nil {
		panic(err)

	}
	session.SetMode(mgo.Monotonic, true)

	defer session.Close()
	session.SetMode(mgo.Monotonic, true)

	file_list := getFilelist(*start_path)
	wg.Add(len(file_list))
	for _, file := range file_list {
		go Md5_Check(&wg, file)

	}
	wg.Wait()
	//	fmt.Println(len(file_list))

}
