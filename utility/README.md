# utility code usage



### requirment

```  
pip install pandas
```

### PepLinkPep2SMILES

peptide 서열과 linker를 연결해서 SMILES로 return 하는 코드
  
- Peptide-linker-Peptide    
- Peptide-linker  
- linker-Peptide  
  
세가지 형태의 연결 지원  

몇가지 검증필요 (proline에 링커 연결 등)
Peptide C-term 말단 Carboxyl 기의 OH 잔기를 NH2로 변환하는 옵션 지원  

linker의 amine기 부분에 아미노산 체인 연결을 기본 전제로 함  
앞부분의 N-term peptide 의 carboxyl기와 linker 연결은 [\*1]로 표현, C-term peptide의 amine기와 linker 연결은 [\*2] 로 표현하여 링커에 결합 위치를 명시해야 함 (example.csv 파일 참고)  
결합 위치를 표현하지 않은경우, linker의 맨앞 혹은 맨뒤 구조 표현식에 결합함

```
python utility/PepLinkPep2SMILES.py help  
python utility/PepLinkPep2SMILES.py -i example_PepLinkPep2SMILES.csv
```
