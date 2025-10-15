#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
링커의 원래 구조를 보존하는 아미노산-링커 구조체 SMILES 변환기
"""

import re,argparse, sys
import pandas as pd

class PeptideLinkerBuilder:
    def __init__(self):
        # 표준 아미노산의 side chain SMILES
        self.amino_acid_sidechains = {
            'A': 'C',                    # Alanine
            'R': 'CCCNC(=N)N',          # Arginine
            'N': 'CC(=O)N',             # Asparagine
            'D': 'CC(=O)O',             # Aspartic acid
            'C': 'CS',                  # Cysteine
            'E': 'CCC(=O)O',            # Glutamic acid
            'Q': 'CCC(=O)N',            # Glutamine
            'G': '',                    # Glycine
            'H': 'CC1=CNC=N1',          # Histidine
            'I': 'C(C)CC',              # Isoleucine
            'L': 'CC(C)C',              # Leucine
            'K': 'CCCCN',               # Lysine
            'M': 'CCSC',                # Methionine
            'F': 'CC1=CC=CC=C1',        # Phenylalanine
            'P': '',                    # Proline - 특별처리
            'S': 'CO',                  # Serine
            'T': 'C(C)O',               # Threonine
            'W': 'CC1=CNC2=CC=CC=C21',  # Tryptophan
            'Y': 'CC1=CC=C(C=C1)O',     # Tyrosine
            'V': 'C(C)C'               # Valine
        }
    
    def validate_sequence(self, sequence):
        """아미노산 서열이 유효한지 확인"""
        # 빈 문자열은 허용
        if not sequence:
            return True
    
        invalid_residues = [aa for aa in sequence if aa not in self.amino_acid_sidechains]
        if invalid_residues:
            raise ValueError(f"유효하지 않은 아미노산: {invalid_residues}")
        return True
    
    def parse_linker_attachment_points(self, linker_smiles):
        """
        링커 SMILES에서 연결점을 파싱
        연결점은 [*1], [*2] 형태로 표시되어야 함
        [*1]은 N-말단, [*2]는 C-말단 연결점
        """
        # [*1]과 [*2]를 각각 찾기
        has_star1 = '[*1]' in linker_smiles
        has_star2 = '[*2]' in linker_smiles
    
        # 연결점이 하나라도 있으면 명시적 연결점으로 처리
        if has_star1 or has_star2:
            return {
                'linker': linker_smiles,
                'n_attachment': '[*1]' if has_star1 else None,
                'c_attachment': '[*2]' if has_star2 else None,
                'has_explicit_points': True
            }
        else:
            return {
                'linker': linker_smiles,
                'n_attachment': None,
                'c_attachment': None,
                'has_explicit_points': False
            }
    
    def build_peptide_chain_for_connection(self, sequence, connection_type):
        """
        연결용 펩타이드 체인 생성
        connection_type: 'n_terminal' (N-말단용), 'c_terminal' (C-말단용)
        """
        if not sequence:
            return ""
        
        if len(sequence) == 1:
            aa = sequence[0]
            sidechain = self.amino_acid_sidechains[aa]
            
            if connection_type == 'n_terminal':
                # N-말단용: 자유 아민으로 시작, 카르보닐 탄소로 끝남
                if aa == 'P':
                    return "N1CCC[CH]1C(=O)"
                elif aa == 'G':
                    return "NCC(=O)"
                else:
                    return f"N[CH]({sidechain})C(=O)"
            else:  # c_terminal
                # C-말단용: 질소로 시작, 자유 카르복실로 끝남
                if aa == 'P':
                    return "N1CCC[CH]1C(=O)O"
                elif aa == 'G':
                    return "NCC(=O)O"
                else:
                    return f"N[CH]({sidechain})C(=O)O"
        
        # 다중 아미노산 체인
        peptide_parts = []
        
        for i, aa in enumerate(sequence):
            sidechain = self.amino_acid_sidechains[aa]
            
            if i == 0:
                # 첫 번째 아미노산
                if connection_type == 'n_terminal':
                    # N-말단용: 자유 아민으로 시작
                    if aa == 'P':
                        peptide_parts.append("N1CCC[CH]1C(=O)")
                    elif aa == 'G':
                        peptide_parts.append("NCC(=O)")
                    else:
                        peptide_parts.append(f"N[CH]({sidechain})C(=O)")
                else:  # c_terminal
                    # C-말단용: 질소로 시작
                    if aa == 'P':
                        peptide_parts.append("N1CCC[CH]1C(=O)")
                    elif aa == 'G':
                        peptide_parts.append("NCC(=O)")
                    else:
                        peptide_parts.append(f"N[CH]({sidechain})C(=O)")
            elif i == len(sequence) - 1:
                # 마지막 아미노산
                if connection_type == 'n_terminal':
                    # N-말단용: 카르보닐 탄소로 끝남
                    if aa == 'P':
                        peptide_parts.append("N1CCC[CH]1C(=O)")
                    elif aa == 'G':
                        peptide_parts.append("NCC(=O)")
                    else:
                        peptide_parts.append(f"N[CH]({sidechain})C(=O)")
                else:  # c_terminal
                    # C-말단용: 자유 카르복실로 끝남
                    if aa == 'P':
                        peptide_parts.append("N1CCC[CH]1C(=O)O")
                    elif aa == 'G':
                        peptide_parts.append("NCC(=O)O")
                    else:
                        peptide_parts.append(f"N[CH]({sidechain})C(=O)O")
            else:
                # 중간 아미노산
                if aa == 'P':
                    peptide_parts.append("N1CCC[CH]1C(=O)")
                elif aa == 'G':
                    peptide_parts.append("NCC(=O)")
                else:
                    peptide_parts.append(f"N[CH]({sidechain})C(=O)")
        
        return ''.join(peptide_parts)

    def convert_last_carboxyl_to_amide(self, smiles):
        """마지막 C(=O)O를 C(=O)N으로 바꾸기"""
        pattern = 'C(=O)O'
        last_index = smiles.rfind(pattern)

        if last_index != -1:
            # C(=O)O의 마지막 O만 N으로 바꾸기
            o_pos = last_index + len(pattern) - 1
            return smiles[:o_pos] + 'N' + smiles[o_pos + 1:]
        return smiles
    
    def connect_with_explicit_points(self, n_peptide_seq, linker_info, c_peptide_seq, remove_n):
	    """
	    명시적 연결점을 사용하여 구조 연결 - 원래 링커 구조 보존
	    링커의 N과 펩타이드의 N 중복 처리
	    remove_n: True이면 C-말단 펩타이드의 첫 N을 제거, False이면 유지
	    
	    지원하는 경우들:
	    1. N-term + C-term 모두 있는 경우
	    2. N-term만 있는 경우 (c_peptide_seq = "")
	    3. C-term만 있는 경우 (n_peptide_seq = "")
	    """
	    linker_smiles = linker_info['linker']
	    n_point = linker_info['n_attachment']
	    c_point = linker_info['c_attachment']
	    
	    # 펩타이드 체인 생성 (빈 서열이면 빈 문자열 반환)
	    n_peptide = self.build_peptide_chain_for_connection(n_peptide_seq, 'n_terminal') if n_peptide_seq else ""
	    c_peptide = self.build_peptide_chain_for_connection(c_peptide_seq, 'c_terminal') if c_peptide_seq else ""
	    
	    # 링커의 연결점을 펩타이드로 치환
	    connected_smiles = linker_smiles
	    
	    # N-말단 연결점 처리
	    if n_point and n_peptide:
	        connected_smiles = connected_smiles.replace(n_point, n_peptide)
	    elif n_point and not n_peptide:
	        # n_peptide가 없으면 [*1]을 그냥 제거
	        connected_smiles = connected_smiles.replace(n_point, '')
	    
	    # C-말단 연결점 처리
	    if c_point and c_peptide:
	        # N[*2] 패턴 처리
	        if remove_n == "True" and 'N' + c_point in linker_smiles:
	            if c_peptide.startswith('N'):
	                c_peptide_without_n = c_peptide[1:]
	                connected_smiles = connected_smiles.replace(c_point, c_peptide_without_n)
	            else:
	                connected_smiles = connected_smiles.replace(c_point, c_peptide)
	        else:
	            connected_smiles = connected_smiles.replace(c_point, c_peptide)
	    elif c_point and not c_peptide:
	        # c_peptide가 없으면 [*2]를 그냥 제거
	        connected_smiles = connected_smiles.replace(c_point, '')
	    
	    return connected_smiles
    
    def build_with_preserved_linker(self, n_terminal_seq, linker_with_points, c_terminal_seq, remove_n):
        """
        링커 구조를 보존하면서 펩타이드 연결
        """
        # 빈 문자열도 허용하도록 검증
        if n_terminal_seq:
            self.validate_sequence(n_terminal_seq)
        if c_terminal_seq:
            self.validate_sequence(c_terminal_seq)
    
            # 링커 연결점 파싱
        linker_info = self.parse_linker_attachment_points(linker_with_points)
    
        if linker_info['has_explicit_points']:
            # 명시적 연결점이 있는 경우 - 원래 구조 보존
            full_smiles = self.connect_with_explicit_points(n_terminal_seq, linker_info, c_terminal_seq, remove_n)
        else:
            # 연결점이 없는 경우 - 기본적인 직접 연결
            full_smiles = self.build_direct_connection(n_terminal_seq, linker_with_points, c_terminal_seq)
    
        return full_smiles
    
    def build_direct_connection(self, n_terminal_seq, linker_smiles, c_terminal_seq):
        """
        연결점 표시가 없는 경우의 직접 연결 (원래 링커 구조 보존)
        """
        n_peptide = self.build_peptide_chain_for_connection(n_terminal_seq, 'n_terminal')
        c_peptide = self.build_peptide_chain_for_connection(c_terminal_seq, 'c_terminal')
        
        # 링커를 중간에 그대로 삽입
        return f"{n_peptide}{linker_smiles}{c_peptide}"
    
    def build_structure_info(self, n_terminal_seq, linker_smiles, c_terminal_seq, NH2_endpoint, remove_n, method="auto"):
        """구조 정보와 함께 결과 반환"""
        try:
            if method == "auto":
                # 연결점 표시자가 있으면 explicit, 없으면 direct
                if '[*' in linker_smiles:
                    smiles = self.build_with_preserved_linker(n_terminal_seq, linker_smiles, c_terminal_seq, remove_n)
                    used_method = "explicit_attachment_preserved"
                else:
                    smiles = self.build_direct_connection(n_terminal_seq, linker_smiles, c_terminal_seq)
                    used_method = "direct_connection_preserved"
            elif method == "preserved":
                smiles = self.build_with_preserved_linker(n_terminal_seq, linker_smiles, c_terminal_seq, remove_n)
                used_method = "preserved_linker_structure"
            else:
                smiles = self.build_direct_connection(n_terminal_seq, linker_smiles, c_terminal_seq)
                used_method = "direct_connection_preserved"
            
            # NH2 endpoint 처리
            if NH2_endpoint == "True":
                smiles = self.convert_last_carboxyl_to_amide(smiles)

            info = {
                'n_terminal_sequence': n_terminal_seq,
                'linker_smiles': linker_smiles,
                'c_terminal_sequence': c_terminal_seq,
                'total_structure_smiles': smiles,
                'method_used': used_method,
                'total_aa_length': len(n_terminal_seq) + len(c_terminal_seq)
            }
            
            return info
            
        except Exception as e:
            return {'error': str(e)}



def main():
    """메인 실행 함수"""
    builder = PeptideLinkerBuilder()



    parser = argparse.ArgumentParser(description = "Peptide-Linker-Peptide name to smiles, see example.csv")
    parser.add_argument("-i","--input",help="example.csv file, do not change header Nterm, Linker, Cterm", default="example.csv")
    parser.add_argument("-E","--Endpoint",help="NH2_endpoint: peptide Cterm Carboxyl O=C-OH to O=C-NH2, True or False", default="True")
    parser.add_argument("-r","--remove_N",help="Cterm peptide amine site remove when amine linking to linker[*2], True or False", default="True")
    
    args = parser.parse_args()

    infile = args.input
    remove_n = args.remove_N
    NH2_endpoint = args.Endpoint


    InputTable = pd.read_csv(infile)

    lists = pd.read_csv(infile, keep_default_na=False).fillna('')

    headers = lists.columns.tolist()
    required = ['Nterm', 'Cterm', 'Linker']
    missing = [col for col in required if col not in headers]
    if missing:
        sys.exit(f"input file error, check input CSV file - {', '.join(missing)}")

    lists['Result'] = ''
    for i in range(lists.__len__()):
        n_term = lists["Nterm"][i]
        linker_explicit = lists["Linker"][i]
        c_term = lists["Cterm"][i]
        #print(n_term)
        #print(linker_explicit)
        #print(c_term)

        result = builder.build_structure_info(n_term, linker_explicit, c_term, NH2_endpoint, remove_n, method="auto")
    
        if 'error' in result:
            print(f"에러: {result['error']}")
        else:
            #print(f"N-말단: {result['n_terminal_sequence']}")
            #print(f"링커: {result['linker_smiles']}")
            #print(f"C-말단: {result['c_terminal_sequence']}")
            #print(f"방법: {result['method_used']}")
            #print(f"결과: {result['total_structure_smiles']}")

            lists.loc[i, 'Result'] = result['total_structure_smiles']
            
    lists.to_csv("result.csv", index=False)
    



if __name__ == "__main__":
    main()