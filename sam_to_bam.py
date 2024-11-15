import os
import subprocess

def convert_sam_to_bam(input_dir, output_dir):
    # 確保輸出目錄存在
    os.makedirs(output_dir, exist_ok=True)
    
    # 獲取所有 .sam 文件
    sam_files = [f for f in os.listdir(input_dir) if f.endswith('.sam')]
    
    for sam_file in sam_files:
        # 構建輸入和輸出文件路徑
        sam_path = os.path.join(input_dir, sam_file)
        bam_file = sam_file.replace('.sam', '.bam')
        bam_path = os.path.join(output_dir, bam_file)
        
        # 構建 samtools 命令
        samtools_command = [
            'samtools', 'view', '-Sb', sam_path, '-o', bam_path
        ]
        
        # 執行命令
        try:
            subprocess.run(samtools_command, check=True)
            print(f"Converted {sam_file} to {bam_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error converting {sam_file}: {e}")

if __name__ == "__main__":
    input_dir = input("Enter the path to the directory containing SAM files: ")
    output_dir = input("Enter the path to the directory where BAM files will be saved: ")
    
    convert_sam_to_bam(input_dir, output_dir)
