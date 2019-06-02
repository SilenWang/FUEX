package main

import (
	"bufio"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"io/ioutil"
	"log"
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/brentp/vcfgo"
	"github.com/urfave/cli"
)

type geneRecord struct {
	exon        [][]int64
	intron      [][]int64
	transRegion [2]int64
	cdsRegion   [2]int64
	chr         string
	strand      string
	geneSymbol  string
}

func reverSlice(a [][]int64) [][]int64 {
	for i := len(a)/2 - 1; i >= 0; i-- {
		opp := len(a) - 1 - i
		a[i], a[opp] = a[opp], a[i]
	}
	return a
}

// 将string array 转数字
func string2int(sa []string) ([]int64, error) {
	si := make([]int64, 0, len(sa))
	for _, a := range sa {
		i, err := strconv.ParseInt(a, 10, 64)
		if err != nil {
			return si, err
		}
		si = append(si, i)
	}
	return si, nil
}

func getDB(dbFile string) (map[string]geneRecord, map[string][]string) {
	tran2info := make(map[string]geneRecord) // map 套 结构体
	gene2tran := make(map[string][]string)   // map 套 string列表
	file, err := os.Open(dbFile)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	buf := bufio.NewReader(file)
	for { // 循环读行
		line, err := buf.ReadString('\n')
		if err != nil || err == io.EOF {
			break
		}
		line = strings.TrimSpace(line)
		cols := strings.Split(line, "\t")

		// 填充 gene2tran
		trans := cols[1]
		gene := cols[12]
		if _, exist := tran2info[gene]; exist { // 这个判断怎么看上去是反了?
			gene2tran[gene] = []string{trans}
		} else {
			gene2tran[gene] = append(gene2tran[gene], trans)
		}

		// 准备填充 tran2info 的信息
		strand := cols[3]
		chr := cols[2]
		start, _ := strconv.ParseInt(cols[4], 10, 64)
		end, _ := strconv.ParseInt(cols[5], 10, 64)
		transRegion := [2]int64{start, end}

		start, _ = strconv.ParseInt(cols[6], 10, 64)
		end, _ = strconv.ParseInt(cols[7], 10, 64)
		cdsRegion := [2]int64{start, end}

		// 解析内含子及外显子
		exonStarts, _ := string2int(strings.Split(cols[9], ","))
		exonEnds, _ := string2int(strings.Split(cols[10], ","))
		exonFrames, _ := string2int(strings.Split(cols[15], ","))

		// fmt.Println(exons)
		exons := [][]int64{}   // Slice 嵌套 Array
		introns := [][]int64{} // Slice 嵌套 Array
		if cols[3] == "+" {
			for idx, _ := range exonStarts {
				exons = append(exons, []int64{exonStarts[idx], exonEnds[idx], exonFrames[idx]})
				if idx != len(exonStarts)-1 {
					introns = append(introns, []int64{exonEnds[idx], exonStarts[idx+1]})
				}
			}
		} else {
			for idx, _ := range exonStarts {
				exons = reverSlice(append(exons, []int64{exonStarts[idx], exonEnds[idx], exonFrames[idx]}))
				if idx != len(exonStarts)-1 {
					introns = reverSlice(append(introns, []int64{exonEnds[idx], exonStarts[idx+1]}))
				}
			}
		}
		// geneRecord填充
		record := geneRecord{
			exon:        exons,
			intron:      introns,
			transRegion: transRegion,
			cdsRegion:   cdsRegion,
			chr:         chr,
			strand:      strand,
			geneSymbol:  gene,
		}
		tran2info[trans] = record
	}
	return tran2info, gene2tran
}

func breakSearch(breakP string, info geneRecord) (string, bool) {
	if chr := strings.Split(breakP, ":")[0]; chr != info.chr {
		return "", false
	}
	pos, _ := strconv.ParseInt(strings.Split(breakP, ":")[1], 10, 64)
	for idx, edges := range info.exon {
		if edges[0] <= pos && pos < edges[1] {
			idxStr := strconv.Itoa(idx + 1)
			return strings.Join([]string{"Ex", idxStr}, ""), true
		}
	}
	for idx, edges := range info.intron {
		if edges[0] <= pos && pos < edges[1] {
			idxStr := strconv.Itoa(idx + 1)
			return strings.Join([]string{"In", idxStr}, ""), true
		}
	}
	return "", false
}

func breakCheck(breakA string, breakB string,
	transA string, transB string,
	tran2info map[string]geneRecord) (map[string]string, string, error) {

	// 在breakA能在transB里找到, 或breakB能在 transA找到时, 两个断点颠倒, 但可能不是目标fusion
	// 在breakA能在transB里找到, 且breakB能在 transA找到时, 两个断点颠倒, 并且是目标fusion
	infoA := tran2info[transA]
	infoB := tran2info[transB]
	locAinA, AinA := breakSearch(breakA, infoA)
	locAinB, AinB := breakSearch(breakA, infoB)
	locBinB, BinB := breakSearch(breakB, infoB)
	locBinA, BinA := breakSearch(breakB, infoA)

	// fmt.Println(locAinA, locAinB, locBinB, locBinA)

	locates := make(map[string]string)

	switch {
	case AinA && BinB: // 目标fusion, 对应正确
		locates["breakA"] = locAinA
		locates["breakB"] = locBinB
		return locates, "right", nil
	case AinB && BinA: // 目标fusion, 对应颠倒
		locates["breakA"] = locAinB
		locates["breakB"] = locBinA
		return locates, "wrong", nil
	default:
		err := errors.New("Not Fusion")
		return locates, "", err
	}
}

func reverse(ori string) string {
	switch {
	case ori == ">":
		ori = "<"
	case ori == "<":
		ori = ">"
	}
	return ori
}

func geneRep(str string, geneA string, geneB string) string {
	str = strings.ReplaceAll(str, "GENE_A", geneA)
	str = strings.ReplaceAll(str, "GENE_B", geneB)
	return str
}

func recordParser(variant *vcfgo.Variant,
	svtype string,
	tran2info map[string]geneRecord,
	gene2tran map[string][]string,
	targetTrans map[string]interface{}) (bndStr string, err error) { //ann 解析的部分直接写这里了
	// 本程序的解析器只从ANN内读取基因注释信息, 因此相较源程序更为通用, 可以将不同类型的解读写到一起

	// 获取类型信息
	// svtype, _ := variant.Info().Get("SVTYPE")

	// 断点及融合形式获取
	var breakA string
	var breakB string
	orientation := [4]string{}
	switch {
	case svtype == "BND":
		breakA = strings.Join([]string{variant.Chromosome, strconv.FormatUint(variant.Pos, 10)}, ":")
		reg := regexp.MustCompile(`CHR.+:\d+`)
		breakB = strings.Replace(reg.FindString(variant.Alt()[0]), "CHR", "chr", 1)
		alt := variant.Alt()[0] // 只会有一个
		switch {
		case strings.Contains(alt, "["):
			alts := strings.Split(alt, "[")
			switch {
			case alts[0] == "":
				// "GENE_B<GENE_A>"
				orientation[0] = "GENE_B"
				orientation[1] = "<"
				orientation[2] = "GENE_A"
				orientation[3] = ">"
			default:
				// "GENE_A>GENE_B>"
				orientation[0] = "GENE_A"
				orientation[1] = ">"
				orientation[2] = "GENE_B"
				orientation[3] = ">"
			}
		case strings.Contains(alt, "]"):
			alts := strings.Split(alt, "]")
			switch {
			case alts[0] == "":
				// "GENE_B>GENE_A>"
				orientation[0] = "GENE_B"
				orientation[1] = ">"
				orientation[2] = "GENE_A"
				orientation[3] = ">"
			default:
				// "GENE_A>GENE_B<"
				orientation[0] = "GENE_A"
				orientation[1] = ">"
				orientation[2] = "GENE_B"
				orientation[3] = "<"
			}
		}
	default:
		posA := variant.Pos
		posB := uint64(variant.End())
		breakA = strings.Join([]string{variant.Chromosome, strconv.FormatUint(posA, 10)}, ":")
		breakB = strings.Join([]string{variant.Chromosome, strconv.FormatUint(uint64(variant.End()), 10)}, ":")
		switch {
		case posA > posB:
			switch {
			case svtype == "INV_1":
				// 修正svtype的值, 防止影响后面输出
				// {GENE_A}>{GENE_B}<
				svtype = "INV"
				orientation[0] = "GENE_A"
				orientation[1] = ">"
				orientation[2] = "GENE_B"
				orientation[3] = "<"
			case svtype == "INV_2":
				// 修正svtype的值, 防止影响后面输出
				// "{GENE_A}<{GENE_B}>"
				svtype = "INV"
				orientation[0] = "GENE_A"
				orientation[1] = "<"
				orientation[2] = "GENE_B"
				orientation[3] = ">"
			case svtype == "DEL":
				// "{GENE_B}>{GENE_A}>"
				orientation[0] = "GENE_B"
				orientation[1] = ">"
				orientation[2] = "GENE_A"
				orientation[3] = ">"
			case svtype == "DUP":
				// "{GENE_A}>{GENE_B}>"
				orientation[0] = "GENE_A"
				orientation[1] = ">"
				orientation[2] = "GENE_B"
				orientation[3] = ">"
			}
		default:
			switch {
			case svtype == "INV_1":
				// 修正svtype的值, 防止影响后面输出
				// "{GENE_B}>{GENE_A}<"
				svtype = "INV"
				orientation[0] = "GENE_B"
				orientation[1] = ">"
				orientation[2] = "GENE_A"
				orientation[3] = "<"
			case svtype == "INV_2":
				// 修正svtype的值, 防止影响后面输出
				// "{GENE_B}<{GENE_A}>"
				svtype = "INV"
				orientation[0] = "GENE_B"
				orientation[1] = "<"
				orientation[2] = "GENE_A"
				orientation[3] = ">"
			case svtype == "DEL":
				// "{GENE_A}>{GENE_B}>"
				orientation[0] = "GENE_A"
				orientation[1] = ">"
				orientation[2] = "GENE_B"
				orientation[3] = ">"
			case svtype == "DUP":
				// "{GENE_B}>{GENE_A}>"
				orientation[0] = "GENE_B"
				orientation[1] = ">"
				orientation[2] = "GENE_A"
				orientation[3] = ">"
			}
		}
	}
	// 获取不需要解析的信息
	pe := variant.Samples[0].Fields["PE"]
	sr := variant.Samples[0].Fields["SR"]

	// 定义准备获取的信息
	var transA, transB string
	var geneA, geneB string

	callID := variant.Id()
	anns, _ := variant.Info().Get("ANN") // 这里返回的是interface, 需要自行把值取出来
	// 这里value可能是string或者[]string, 同时只能取其中一种, 所以加入判断
	var annTarget string
	if annObj, ok := anns.(string); ok { // 单个注释, 直接解析
		annTarget = annObj
	} else { // 多个注释, 暂时取第一个, 因为这份程序里通过基因解析转录本, 所以只需要第一个就行(除非对应了多个基因?)
		annObj, _ := anns.([]string)
		annTarget = annObj[0]
	}
	annCell := strings.Split(annTarget, "|")
	if strings.Contains(annCell[3], "&") {
		genes := strings.Split(annCell[3], "&")
		// 获取gene symbol
		geneA = genes[0]
		geneB = genes[1]

		// 获取转录本ID
		if value := targetTrans[geneA]; value != nil {
			transA = targetTrans[geneA].(string)
		} else {
			if value := gene2tran[geneA]; value != nil { // 可能存在gene symbol在文件库中不存在, 不再解析
				transA = gene2tran[geneA][0] // 直接取第一个
			} else {
				fmt.Printf("WARING: Gene Symbol not in database, caller id: %s \n", callID)
				err = errors.New("No Gene Symbol Found")
				return "", err
			}
		}
		if value := targetTrans[geneB]; value != nil {
			transB = targetTrans[geneB].(string)
		} else {
			if value := gene2tran[geneB]; value != nil { // 可能存在gene symbol在文件库中不存在, 不再解析
				transB = gene2tran[geneB][0] // 直接取第一个
			} else {
				fmt.Printf("WARING: Gene Symbol not in database, caller id: %s \n", callID)
				err = errors.New("No Gene Symbol Found")
				return "", err
			}
		}
	} else {
		fmt.Println("WARING: Fusion within single gene in %s, passed", callID)
		err = errors.New("Fusion within single gene")
		return "", err
	}

	// 获取转录方向
	strandA := tran2info[transA].strand
	strandB := tran2info[transB].strand

	// 到目前位置gene/trans和 break是没有对应的, 现在进行修正
	// AB以break情况为准, 所以调换的是gene和trans
	locates, consist, err := breakCheck(breakA, breakB, transA, transB, tran2info)
	// 由于map不能方便的深复制, 如果需要原始信息, 需要在此生成

	if err != nil {
		err := errors.New("Not Fusion")
		return "", err
	}
	if consist == "wrong" {
		geneA, geneB = geneB, geneA
		transA, transB = transB, transA
		strandA, strandB = strandB, strandA
	}

	// 基本信息准备完毕,  进行方向矫正
	modOri := orientation
	if strandA == "-" {
		for idx, val := range modOri {
			if val == "GENE_A" {
				modOri[idx+1] = reverse(modOri[idx+1])
			}
		}
	}
	if strandB == "-" {
		for idx, val := range modOri {
			if val == "GENE_B" {
				modOri[idx+1] = reverse(modOri[idx+1])
			}
		}
	}

	// 判断条目是否真fusion
	var fusion bool
	if modOri[1] == modOri[3] {
		fusion = true
	} else {
		fusion = false
	}

	if fusion {
		// locates外显子号码矫正
		for idx, val := range modOri {
			if strings.Contains(val, "GENE") {
				key := strings.Join([]string{"break", val[len(val)-1:]}, "") // 拼接出key
				locate := locates[key]
				switch {
				case strings.Contains(locate, "In"): // 内含子则校正为外显子
					num, _ := strconv.Atoi(locate[2:])
					if modOri[idx+1] == ">" { // <的话外显子号等于内含子号
						num = num + 1
					}
					// switch {
					// case modOri[idx+1] == ">":
					// 	num = num + 1
					// case modOri[idx+1] == "<":
					// 	num = num - 1
					// }
					// 修正结果替换
					locates[key] = strings.Join([]string{"Ex", strconv.Itoa(num)}, "")
				case strings.Contains(locate, "Ex"): // 断点在Ex则加*提示
					locates[key] = strings.Join([]string{locate, "*"}, "")
				}
			}
		}

		// 所有信息都齐了, 开始准备输出用信息

		var fuseGene string
		switch {
		case modOri[1] == ">":
			fuseGene = strings.Join([]string{modOri[0], modOri[2]}, "-")
		case modOri[1] == "<":
			fuseGene = strings.Join([]string{modOri[2], modOri[0]}, "-")
		}
		fuseGene = geneRep(fuseGene, geneA, geneB)

		var fuseDetail string
		transMap := map[string]string{
			"GENE_A": transA,
			"GENE_B": transB,
		}
		locMap := map[string]string{
			"GENE_A": locates["breakA"],
			"GENE_B": locates["breakB"],
		}
		switch {
		case modOri[1] == ">":
			fuseDetail = strings.Join([]string{
				strings.Join([]string{modOri[0], transMap[modOri[0]], locMap[modOri[0]]}, ":"),
				strings.Join([]string{modOri[2], transMap[modOri[2]], locMap[modOri[2]]}, ":"),
			}, ">")
		case modOri[1] == "<":
			fuseDetail = strings.Join([]string{
				strings.Join([]string{modOri[2], transMap[modOri[2]], locMap[modOri[2]]}, ":"),
				strings.Join([]string{modOri[0], transMap[modOri[0]], locMap[modOri[0]]}, ":"),
			}, ">")
		}
		fuseDetail = geneRep(fuseDetail, geneA, geneB)

		outStrList := []string{
			fuseGene,                            // fusion
			fuseDetail,                          // fusion_detail
			svtype,                              // svtype
			callID,                              // caller_id
			breakA,                              // break_A
			breakB,                              // break_B
			strings.Join([]string{pe, sr}, ":"), // PE:SR
		}

		outLine := strings.Join(outStrList, "\t")
		// fmt.Println(outLine)
		// fmt.Println(fuseGene, fuseDetail, breakA, breakB, strandA, strandB, svtype, pe, sr, locates, orientation, modOri, fusion)
		return outLine, nil
	}
	err = errors.New("Not Invalid fusion")
	return "", err
}

func run(vcfFile string, outFile string) {
	tran2info, gene2tran := getDB("refGene.txt") // 获取数据结构
	// 载入目标转录本信息
	tarTranFile := "RXA_gene_trans.json"
	data, err := ioutil.ReadFile(tarTranFile)
	if err != nil {
		return
	}
	var targetTransFile interface{}
	err = json.Unmarshal(data, &targetTransFile)
	targetTrans := targetTransFile.(map[string]interface{})

	// 准备写入对象
	out, err := os.Create(outFile)
	if err != nil {
		log.Fatal(err)
	}
	out.Close()
	out, err = os.OpenFile(outFile, os.O_WRONLY, 0666)
	if err != nil {
		log.Fatal(err)
	}
	defer out.Close()
	// 为这个文件创建buffered writer
	bufferedWriter := bufio.NewWriter(out)

	// const vcfFile = "test.vcf"
	// 打开文件获取阅读对象
	file, err := os.Open(vcfFile)
	defer file.Close() // 自动关闭
	vcfHandle := io.Reader(file)
	vcfReader, err := vcfgo.NewReader(vcfHandle, false)
	if err != nil {
		panic(err)
	}
	for { // 永循环, 读不到的时候退出
		variant := vcfReader.Read()
		if variant == nil {
			break
		} else {
			sec, _ := variant.Info().Get("SECONDARY")
			if sec == true { // 这里返回类型是interface, 必须判断, 另外如果是sec不进行解析
				break
			}
			svtype, _ := variant.Info().Get("SVTYPE")
			switch {
			case svtype == "INV": // INV两端的端点是两种融合形式, 所以这里特殊处理
				outLine, fuse := recordParser(variant, "INV_1", tran2info, gene2tran, targetTrans)
				if fuse == nil {
					_, err := bufferedWriter.WriteString(
						strings.Join([]string{outLine, "\n"}, ""),
					)
					if err != nil {
						log.Fatal(err)
					}
					bufferedWriter.Flush()
				}
				outLine, fuse = recordParser(variant, "INV_2", tran2info, gene2tran, targetTrans)
				if fuse == nil {
					_, err := bufferedWriter.WriteString(
						strings.Join([]string{outLine, "\n"}, ""),
					)
					if err != nil {
						log.Fatal(err)
					}
					bufferedWriter.Flush()
				}
			// println(outLine)				continue
			default:
				outLine, fuse := recordParser(variant, svtype.(string), tran2info, gene2tran, targetTrans)
				if fuse == nil {
					_, err := bufferedWriter.WriteString(
						strings.Join([]string{outLine, "\n"}, ""),
					)
					if err != nil {
						log.Fatal(err)
					}
					bufferedWriter.Flush()
				}
				// println(outLine)
			}
		}
	}
	vcfReader.Clear()
}

func argParse() {
	// 定义获取的变量
	var vcfFile string
	var out string

	// 指定程序基本信息
	app := cli.NewApp()
	app.Name = "FUEX"
	app.Usage = "Fusion extracter of SV detection tools"
	app.Version = "0.1.1"

	app.Flags = []cli.Flag{
		cli.StringFlag{
			Name:        "i",
			Usage:       "Input file of `VCF` format",
			Destination: &vcfFile,
		},
		cli.StringFlag{
			Name:        "o",
			Usage:       "`OUTPUT` file",
			Destination: &out,
		},
	}

	app.Action = func(c *cli.Context) error {
		fmt.Println(vcfFile, out)
		run(vcfFile, out)
		return nil
	}

	err := app.Run(os.Args)
	if err != nil {
		log.Fatal(err)
	}
}

func main() {
	argParse()
}
