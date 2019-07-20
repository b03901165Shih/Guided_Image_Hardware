# Guided_Image_Hardware
This is a Verilog implemented hardware architecture for 30fps full-HD gray-scale guided image filter

## Usage

1. Run gen_patt_integer_image.m for test pattern gereration
  - Change stripe_ind (0-15) to test different image stripe
  
2. Copy DW_div_function.inc from your designware directory to this folder

3. Run simulation:

  ncverilog tb_IIE.v Guided_Image.v +access+r -y \<DesignWare directory> +libext+.v
  - \<DesignWare directory>: your designware directory path (e.g. /usr/cad/synopsys/synthesis/cur/dw/sim_ver/)

## Visualize

Run read_verilog_image.m for visulization

## Result
<p align="left">
  <img width="100%" height="100%" src="https://github.com/b03901165Shih/Guided_Image_Hardware/blob/master/result/Our-result.png" />
</p>

## Reference

[1] He, Kaiming, Jian Sun, and Xiaoou Tang. "Guided image filtering." IEEE transactions on pattern analysis and machine intelligence 35.6 (2012): 1397-1409.

[2] Kao, Chieh-Chi, Jui-Hsin Lai, and Shao-Yi Chien. "VLSI architecture design of guided filter for 30 frames/s full-HD video." IEEE Transactions on Circuits and Systems for Video Technology 24.3 (2013): 513-524.

