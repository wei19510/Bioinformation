import os
import subprocess
import argparse

def main(input_dir, output_dir):
    # 設定 Q30 和 P40 標準的修剪參數
    quality_cutoff = 30
    min_length = 40

    # 設定引物序列
    forward_primer = 'GTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    reverse_primer = 'CAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAC'

    # 確保輸出目錄存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f'Created output directory: {output_dir}')

    # 瀏覽目錄中的所有文件
    for file_name in os.listdir(input_dir):
        if file_name.endswith('_R1_001.fastq.gz'):
            base_name = file_name.replace('_R1_001.fastq.gz', '')
            r1_file = os.path.join(input_dir, file_name)
            r2_file = os.path.join(input_dir, file_name.replace('_R1_001.fastq.gz', '_R2_001.fastq.gz'))
            
            if os.path.exists(r2_file):
                # 設定輸出文件名
                output_r1 = os.path.join(output_dir, base_name + '_R1_trimmed.fastq.gz')
                output_r2 = os.path.join(output_dir, base_name + '_R2_trimmed.fastq.gz')

                # 使用 cutadapt 進行修剪和引物去除
                cmd = [
                    'cutadapt',
                    '-q', str(quality_cutoff),  # 設定 Q 值
                    '-m', str(min_length),      # 設定最小讀長
                    '-a', forward_primer,       # 去除 forward 引物
                    '-A', reverse_primer,       # 去除 reverse 引物
                    '-o', output_r1,            # 輸出修剪後的 R1 文件
                    '-p', output_r2,            # 輸出修剪後的 R2 文件
                    r1_file,                    # 輸入 R1 文件
                    r2_file                     # 輸入 R2 文件
                ]
                
                # 執行命令
                subprocess.run(cmd, check=True)

                print(f'Processed: {base_name}')
            else:
                print(f'Warning: {r2_file} not found for {r1_file}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim FASTQ files and remove primers.')
    parser.add_argument('input_dir', type=str, help='Path to the input directory containing FASTQ files.')
    parser.add_argument('output_dir', type=str, help='Path to the output directory for trimmed FASTQ files.')

    args = parser.parse_args()

    main(args.input_dir, args.output_dir)
