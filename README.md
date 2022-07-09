# Bioinformatics Course Coding Assignments

## [Semi-Global Alignment](https://github.com/Amirhossein-Rajabpour/Bioinformatics-Course-Assignments/blob/main/HW2/main.py)
In this assignment `semi-global alignment` algorithm was implemented from scratch. `PAM250` matric was used for scoring. <br>
The program takes two strings as input and calculate the alignment score and two modified stirngs after alignment. <br>

Sample input:
```
HEAGAWGHE
PAWHEA
```

Sample output:
```
20
HEAGAWGHE-
---PAW-HEA
```

## [MSA - Star Alignment](https://github.com/Amirhossein-Rajabpour/Bioinformatics-Course-Assignments/blob/main/HW3/main.py)
This code comprises two parts. In the first part `Star Alignment` which is a `Multuple Sequence Alignment` was implemented. Second part was a modification step which had to find blocks that can be modified and replace these blocks in the sequences after alignments. <br>
It takes a number which indicated how many sequences we have. After that sequences are given. It outputs the final MSA score and aligned and modified sequences. <br>

Sample input:
```
4 
TYIMREAQYESAQ
TCIVMREAYE
YIMQEVQQER
WRYIAMREQYES
```

Sample output:
```
51
-TYI-MREAQYESAQ
-TCIVMREA-YE---
--YI-MQEVQQER--
WRYIAMRE-QYES--
```

## [MSA Profile](https://github.com/Amirhossein-Rajabpour/Bioinformatics-Course-Assignments/blob/main/HW4/main.py)
This assignment is mainly about creating a `profile` for `MSA`. Then finding the best subsequence of a given string using profile scores calculated earlier. <br>
The program takes an interger which shows how many sequences are inside MSA and then takes the sequences. after the sequences it takes a long string to find the best subsequence within it. It outputs the best sbsequence found. <br>

Sample input:
```
4
T-CT
--CT
A-CT
ATCT
ATCCTATATCTTCTCTATACTATCCTTCA
```

Sample output:
```
A-CT
```
