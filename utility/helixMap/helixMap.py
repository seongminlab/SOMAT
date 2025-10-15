import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

def get_amino_acid_properties():
    """아미노산의 화학적 속성 분류"""
    return {
        'polar_basic': ['R', 'K', 'H'],
        'polar_acidic': ['D', 'E'],
        'polar_uncharged': ['S', 'T', 'N', 'Q', 'C', 'Y'],
        'nonpolar': ['A', 'V', 'L', 'I', 'M', 'F', 'W', 'P', 'G']
    }

def get_color(amino_acid, properties=None):
    """아미노산 속성에 따른 색상 반환"""
    if properties is None:
        properties = get_amino_acid_properties()
    
    if amino_acid in properties['polar_basic']:
        return '#4169E1'
    elif amino_acid in properties['polar_acidic']:
        return '#DC143C'
    elif amino_acid in properties['polar_uncharged']:
        return '#32CD32'
    elif amino_acid in properties['nonpolar']:
        return '#FFD700'
    else:
        return '#CCCCCC'

def draw_helical_wheel(sequence, title=None, show_connections=True, 
                       save_path=None, figsize=(12, 12), spiral=True):
    """
    펩타이드 서열의 helical wheel 다이어그램 생성
    
    Parameters:
    -----------
    sequence : str
        아미노산 서열
    title : str, optional
        그림 제목
    show_connections : bool, default=True
        연결선 표시 여부
    save_path : str, optional
        저장할 파일 경로
    figsize : tuple, default=(12, 12)
        그림 크기
    spiral : bool, default=True
        나선형 배치 사용 여부 (긴 서열에 권장)
    """
    angle_per_residue = 100  # degrees
    properties = get_amino_acid_properties()
    
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    
    positions = []
    
    # 나선형 배치를 위한 반지름 증가 설정
    if spiral:
        # 18개마다 한 바퀴 도는 것을 고려하여 반지름 증가
        base_radius = 1.0
        radius_increment = 0.15  # 한 바퀴당 증가량
    else:
        base_radius = 1.0
        radius_increment = 0.0
    
    # 각 아미노산 그리기
    for i, aa in enumerate(sequence):
        # 극좌표 계산
        angle = np.radians(i * angle_per_residue)
        
        # 나선형: 반지름을 점진적으로 증가
        if spiral:
            # 18개마다 한 바퀴 = 3.6 residues/turn
            turns = i / 3.6
            radius = base_radius + (radius_increment * turns)
        else:
            radius = base_radius
        
        # 직교좌표 변환
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        positions.append((x, y))
        
        # 색상 결정
        color = get_color(aa, properties)
        
        # 원 크기 조절 (긴 서열일수록 작게)
        circle_radius = max(0.12, 0.18 - (len(sequence) / 500))
        
        # 원 그리기
        circle = Circle((x, y), circle_radius, color=color, 
                       ec='black', linewidth=2, zorder=3)
        ax.add_patch(circle)
        
        # 아미노산 문자 표시
        fontsize = max(10, 16 - (len(sequence) / 50))
        ax.text(x, y, aa, ha='center', va='center', 
                fontsize=fontsize, fontweight='bold', zorder=4)
        
        # 위치 번호 표시
        offset = circle_radius + 0.15
        label_x = x + offset * np.cos(angle)
        label_y = y + offset * np.sin(angle)
        label_fontsize = max(8, 12 - (len(sequence) / 60))
        ax.text(label_x, label_y, str(i+1), 
                ha='center', va='center', fontsize=label_fontsize, 
                color='gray', zorder=2)
    
    # 연결선 그리기
    if show_connections and len(positions) > 1:
        for i in range(len(positions) - 1):
            x1, y1 = positions[i]
            x2, y2 = positions[i + 1]
            ax.plot([x1, x2], [y1, y2], 'gray', linewidth=0.8, 
                    alpha=0.4, zorder=1)
    
    # 축 범위 자동 조정
    max_radius = base_radius + (radius_increment * len(sequence) / 3.6)
    limit = max_radius + 0.8
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.axis('off')
    
    # 범례 추가
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#4169E1', edgecolor='black', label='Polar / Basic'),
        Patch(facecolor='#DC143C', edgecolor='black', label='Polar / Acidic'),
        Patch(facecolor='#32CD32', edgecolor='black', label='Polar / Uncharged'),
        Patch(facecolor='#FFD700', edgecolor='black', label='Nonpolar'),
        Patch(facecolor='#CCCCCC', edgecolor='black', label='Unknown')
    ]
    ax.legend(handles=legend_elements, loc='upper right', 
              fontsize=11, frameon=True)
    
    # 제목
    if title:
        plt.title(title, fontsize=16, fontweight='bold', pad=20)
    else:
        plt.title(f'Helical Wheel: {sequence[:20]}{"..." if len(sequence) > 20 else ""}', 
                  fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    # 저장
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"그림이 저장되었습니다: {save_path}")
    
    plt.show()

# 사용 예시
if __name__ == "__main__":

    sequence = "VDDALALRPRGPEPAEGLRD"

    if len(sequence) > 18:
        spiral_layout = True
    else:
        spiral_layout = False

    draw_helical_wheel(
        sequence=sequence,
        title=sequence,
        spiral=spiral_layout,
        save_path="helical_wheel_spiral.png"
    )