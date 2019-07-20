module Guided_Image 
#(parameter STRIPEWIDTH = 120, parameter ALPHA = 15, parameter A_FRACBITS = 8, parameter B_FRACBITS = 1)	//S = 2*alpha+1
(
clk, rst_n, Iq, Is, pq, ps, valid_in, ready, read_addr_q, read_addr_s, 
//Second Stage
ak, bk, valid_inter_sigab, 
//Third Stage
aq, bq, Sigak, Sigbk, read_addr_q_ab, 
//Fourth Stage
I_in, qi, read_addr_Iin, valid_out);

//First Stage
parameter SIGI_BITS 	= $clog2((2*ALPHA+1)*(2*ALPHA+1))+8;	//18
parameter SIGII_BITS 	= $clog2((2*ALPHA+1)*(2*ALPHA+1))+16;	//26
//Second Stage
parameter ADDR_BITS 		= $clog2(1080*(STRIPEWIDTH+4*ALPHA));	//18
//Third Stage
parameter ADDR_BITS_THIRD 	= $clog2((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA));	//13
parameter SIGIA_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+(4+A_FRACBITS+1);	//10+13=23
parameter SIGIB_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+(12+B_FRACBITS+1);	//10+14=24

input 				clk;
input 				rst_n;
input 				valid_in;
input [7:0] 		Iq;
input [7:0] 		Is;
input [7:0] 		pq;
input [7:0] 		ps;
input				ready; //indicate testbench acknowledge read_addr

input signed [4+A_FRACBITS:0]	aq;
input signed [12+B_FRACBITS:0]	bq;

input [7:0]			I_in;

//Read Iq/pq, Is/ps
output [ADDR_BITS-1:0]  read_addr_q;
output [ADDR_BITS-1:0]  read_addr_s;
output [ADDR_BITS_THIRD-1:0]	read_addr_q_ab;
output [ADDR_BITS-1:0]  read_addr_Iin;

output signed [4+A_FRACBITS:0]	ak;
output signed [12+B_FRACBITS:0]	bk;
output 					valid_out;
output [7:0]			qi;

output signed [SIGIA_BITS-1:0]	Sigak;
output signed [SIGIB_BITS-1:0]	Sigbk;
output							valid_inter_sigab;

wire [SIGI_BITS-1:0]	SigIi, SigIi_pip;
wire [SIGI_BITS-1:0]	Sigpi, Sigpi_pip;
wire [SIGII_BITS-1:0]	SigIiIi;
wire [SIGII_BITS-1:0]	SigIipi;
wire					valid_inter, valid_inter_bk, valid_inter_q;


wire signed [4+A_FRACBITS:0]	ak_inter;


First_Stage #(STRIPEWIDTH, ALPHA) first_stage_u
(
.clk(clk), .rst_n(rst_n), .Iq(Iq), .Is(Is), .pq(pq), .ps(ps), .valid_in(valid_in), 
.SigIi(SigIi), .Sigpi(Sigpi), .SigIiIi(SigIiIi), .SigIipi(SigIipi), .valid_out(valid_inter),
.ready(ready), .read_addr_q(read_addr_q), .read_addr_s(read_addr_s)
);

Second_Stage_AK #(STRIPEWIDTH, ALPHA, A_FRACBITS)	second_stage_ak_u
(
.clk(clk), .rst_n(rst_n), .SigIi(SigIi), .Sigpi(Sigpi), .SigIiIi(SigIiIi), .SigIipi(SigIipi), .valid_in(valid_inter), 
.ak(ak_inter), .SigIi_pip(SigIi_pip), .Sigpi_pip(Sigpi_pip), .valid_out(valid_inter_bk)
);

Second_Stage_BK #(STRIPEWIDTH, ALPHA, A_FRACBITS, B_FRACBITS)	second_stage_bk_u
(
.clk(clk), .rst_n(rst_n), .ak_in(ak_inter), .SigIi(SigIi_pip), .Sigpi(Sigpi_pip), .valid_in(valid_inter_bk), 
.ak(ak), .bk(bk), .valid_out(valid_inter_sigab)
);

Third_Stage #(STRIPEWIDTH, ALPHA, A_FRACBITS, B_FRACBITS) third_stage_ab
(
.clk(clk), .rst_n(rst_n), .aq(aq), .as(ak), .bq(bq), .bs(bk), .valid_in(valid_inter_sigab), 
.Sigak(Sigak), .Sigbk(Sigbk), .valid_out(valid_inter_q), .read_addr_q(read_addr_q_ab)
);

Fourth_Stage #(STRIPEWIDTH, ALPHA, A_FRACBITS , B_FRACBITS) fourth_stage_u
(
.clk(clk), .rst_n(rst_n), .Sigak(Sigak), .Sigbk(Sigbk), .I_in(I_in), 
.read_addr_Iin(read_addr_Iin), .valid_in(valid_inter_q), .qi(qi), .valid_out(valid_out)
);

endmodule

//===========================================================================//
module Fourth_Stage
#(parameter STRIPEWIDTH = 120, parameter ALPHA = 15, parameter A_FRACBITS = 8, parameter B_FRACBITS = 1)	//S = 2*alpha+1
(clk, rst_n, Sigak, Sigbk, I_in, read_addr_Iin, valid_in, qi, valid_out);

parameter AREA_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1));	//10

//Need to read from (1080)*(STRIPEWIDTH+4*ALPHA) I_in memory
parameter ADDR_BITS = $clog2(1080*(STRIPEWIDTH+4*ALPHA));	//18

parameter SIGIA_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+(4+A_FRACBITS+1);	//10+13=23
parameter SIGIB_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+(12+B_FRACBITS+1);	//10+14=24

parameter NUMOFROWS = 1080;	//1080p :^)
parameter NUMOFCOLS = STRIPEWIDTH;	//ex: 1120

parameter ROWBITS = $clog2(NUMOFROWS);
parameter COLBITS = $clog2(NUMOFCOLS);

parameter NUM_STAGES = 4;

input 				clk;
input 				rst_n;
input 				valid_in;
input signed [SIGIA_BITS-1:0]	Sigak;
input signed [SIGIB_BITS-1:0]	Sigbk;
input [7:0]	I_in;

output valid_out;
output reg [7:0]	qi;
output reg [ADDR_BITS-1:0] read_addr_Iin;

reg [ADDR_BITS-1:0]  read_addr_Iin_nxt;// = (!valid_in)?read_addr_Iin:((read_addr_Iin == ((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA)-1))?0:read_addr_Iin+1);


//wire signed [SIGI_BITS:0]	s_SigIi = {1'b0, SigIi};
wire signed [8:0] s_I_in = {1'b0, I_in};
wire signed [SIGIA_BITS+7:0] Sigbk_shift = {Sigbk, {(A_FRACBITS-B_FRACBITS){1'b0}}};

reg [ROWBITS-1:0] row_count, row_count_nxt;
reg [COLBITS-1:0] col_count, col_count_nxt;

reg [AREA_BITS-1:0] area_count, area_count_nxt;


wire signed [SIGIA_BITS+7:0] IiSigak = Sigak*s_I_in;	//should be SIGIA_BITS+8 ???

wire signed [SIGIA_BITS+7:0] Numer = IiSigak+Sigbk_shift; //should be SIGIA_BITS+8 ???

wire signed [AREA_BITS:0] Denom = {1'b0, area_count};

///////////////////////////////////////////////////////////////////////////////////////////
//

reg signed [SIGIA_BITS+7:0] div_a;
reg signed [AREA_BITS:0] div_b;

wire signed [SIGIA_BITS+7:0] quotient_inst; //Need to remove the last (A_FRACBITS) bits

reg signed  [7:0]	qk_nxt;

reg div_in_en;

always@(*)begin
	row_count_nxt = row_count;
	col_count_nxt = col_count;
	area_count_nxt = area_count;
	//output to 0 if negative; output 255 if over 255 and not negative; 4 left 5 in
	qk_nxt = (quotient_inst[SIGIA_BITS+7])?8'd0:((quotient_inst[A_FRACBITS+8])?8'd255:((quotient_inst[A_FRACBITS-1])?quotient_inst[A_FRACBITS+7:A_FRACBITS]+1:quotient_inst[A_FRACBITS+7:A_FRACBITS]));
	read_addr_Iin_nxt = read_addr_Iin;
	if(valid_in) begin
		row_count_nxt = (col_count==(NUMOFCOLS-1))?((row_count==(NUMOFROWS-1))?0:row_count+1):row_count;
		col_count_nxt = (col_count==(NUMOFCOLS-1))?0:col_count+1;
		area_count_nxt = (col_count==(NUMOFCOLS-1))?((row_count<=(ALPHA-1))?area_count+(2*ALPHA+1):((row_count>=(1079-ALPHA))?area_count-(2*ALPHA+1):area_count)):area_count;
		area_count_nxt = (col_count==(NUMOFCOLS-1) && row_count==(NUMOFROWS-1))?(ALPHA+1)*(2*ALPHA+1):area_count_nxt;
		read_addr_Iin_nxt = (col_count==(NUMOFCOLS-1))?((row_count==(NUMOFROWS-1))?0:read_addr_Iin+1+4*ALPHA):read_addr_Iin+1;
	end
end

always@(posedge clk or negedge rst_n)begin
	if(!rst_n)begin
		row_count <= 0;
		col_count <= 0;
		area_count <= (ALPHA+1)*(2*ALPHA+1);
		div_in_en <= 0;
		qi <= 0;
		div_a <= 0;
		div_b <= 0;
		read_addr_Iin <= 30;
	end	else begin
		row_count <= row_count_nxt;
		col_count <= col_count_nxt;	
		area_count <= area_count_nxt;
		div_in_en <= valid_in;
		qi <= qk_nxt;
		div_a <= Numer;
		div_b <= Denom;
		read_addr_Iin <= read_addr_Iin_nxt;
	end
end

//23+8=31 / 10+1 = 11
DW_div_pipe #(SIGIA_BITS+8, AREA_BITS+1, 1, 1, NUM_STAGES, 1, 1, 0)  U1 
(
.clk(clk),   .rst_n(rst_n),   .en(1'b1),
.a(div_a),   .b(div_b),   .quotient(quotient_inst),
.remainder(),   .divide_by_0() 
);

Shift_Reg #(NUM_STAGES, 1) valid_sr
(.clk(clk), .rst_n(rst_n), .in_data(div_in_en), .out_data(valid_out));

endmodule

//==============================================================================//

//===========================================================================//
module Third_Stage
#(parameter STRIPEWIDTH = 120, parameter ALPHA = 15, parameter A_FRACBITS = 8, parameter B_FRACBITS = 1)	//S = 2*alpha+1
(clk, rst_n, aq, as, bq, bs, valid_in, Sigak, Sigbk, valid_out, read_addr_q);

parameter IA_BITS 	= $clog2((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA))+(4+A_FRACBITS+1);	//13+13=26
parameter IB_BITS 	= $clog2((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA))+(12+B_FRACBITS+1);	//13+14=27

//only need to read (2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA) memory (aq/bq)
parameter ADDR_BITS = $clog2((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA));	//13

parameter SIGIA_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+(4+A_FRACBITS+1);	//10+13=23
parameter SIGIB_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+(12+B_FRACBITS+1);	//10+14=24
parameter WINDOW_S = (2*ALPHA+1);

parameter NUMOFROWS = 1080+ALPHA;			//1080p :^)
parameter NUMOFCOLS = STRIPEWIDTH+2*ALPHA;	//ex: 150

parameter ROWBITS = $clog2(NUMOFROWS);
parameter COLBITS = $clog2(NUMOFCOLS);

input 				clk;
input 				rst_n;
input 				valid_in;

input signed [4+A_FRACBITS:0]	aq;
input signed [12+B_FRACBITS:0]	bq;
input signed [4+A_FRACBITS:0]	as;
input signed [12+B_FRACBITS:0]	bs;

//Read aq/bq
output reg [ADDR_BITS-1:0]  read_addr_q;

wire [ADDR_BITS-1:0]  read_addr_q_nxt = (!valid_in)?read_addr_q:((read_addr_q == ((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA)-1))?0:read_addr_q+1);

output signed [SIGIA_BITS-1:0]	Sigak;
output signed [SIGIB_BITS-1:0]	Sigbk;
output 					valid_out;


wire 				 wen, ren1, ren2;
wire [COLBITS-1:0] 	 waddr, raddr1, raddr2;
wire [IA_BITS-1:0] 	 wdata_ai, rdata1_ai, rdata2_ai;
wire [IB_BITS-1:0] 	 wdata_bi, rdata1_bi, rdata2_bi;

reg valid_in_pip;
reg [4+A_FRACBITS:0]  aq_pip, as_pip;
reg [12+B_FRACBITS:0] bq_pip, bs_pip;

//always@(*)begin	
//end


always@(posedge clk or negedge rst_n)begin
	if(!rst_n)begin
		valid_in_pip <= 0;
		aq_pip <= 0;
		as_pip <= 0;
		bq_pip <= 0;
		bs_pip <= 0;
		read_addr_q <= 0;
	end	else begin
		valid_in_pip <= valid_in;
		aq_pip <= aq;
		as_pip <= as;
		bq_pip <= bq;
		bs_pip <= bs;
		read_addr_q <= read_addr_q_nxt;
	end
end


IIE_signed #((4+A_FRACBITS+1),  ALPHA, IA_BITS, SIGIA_BITS, NUMOFROWS, NUMOFCOLS) IIE_ai
(
.clk(clk), .rst_n(rst_n), .Iq(aq_pip), .Is(as_pip), .valid_in(valid_in_pip), .SigIi(Sigak), .valid_out(valid_out), 
.wen(wen), .ren1(ren1), .ren2(ren2), .waddr(waddr), .raddr1(raddr1), .raddr2(raddr2), .wdata(wdata_ai), .rdata1(rdata1_ai), .rdata2(rdata2_ai)
);

IIE_signed #((12+B_FRACBITS+1), ALPHA, IB_BITS, SIGIB_BITS, NUMOFROWS, NUMOFCOLS) IIE_bi
(
.clk(clk), .rst_n(rst_n), .Iq(bq_pip), .Is(bs_pip), .valid_in(valid_in_pip), .SigIi(Sigbk), .valid_out(), 
.wen(), .ren1(), .ren2(), .waddr(), .raddr1(), .raddr2(), .wdata(wdata_bi), .rdata1(rdata1_bi), .rdata2(rdata2_bi)
);

//26+27=53, words: 150
memory_wrapper  #(IA_BITS+IB_BITS, NUMOFCOLS) rf_mem
(
	.clk(clk), .rst_n(rst_n), 
	.wen(wen), .waddr(waddr), .wdata({wdata_ai, wdata_bi}), 
	.ren1(ren1), .raddr1(raddr1), .rdata1({rdata1_ai, rdata1_bi}), 
	.ren2(ren2), .raddr2(raddr2), .rdata2({rdata2_ai, rdata2_bi})
);

endmodule

//===========================================================================//

//===========================================================================//
module IIE_signed
#(parameter I_BITS = 8, parameter ALPHA = 15, parameter II_BITS = 21, parameter SIGI_BITS = 18, parameter NUMOFROWS =1095, parameter NUMOFCOLS = 180)	//S = 2*alpha+1
(clk, rst_n, Iq, Is, valid_in, SigIi, valid_out, wen, ren1, ren2, waddr, raddr1, raddr2, wdata, rdata1, rdata2);

//parameter II_BITS 	= $clog2((2*ALPHA+1)*(STRIPEWIDTH+4*ALPHA))+I_BITS;	//21
//parameter SIGI_BITS	= $clog2((2*ALPHA+1)*(2*ALPHA+1))+I_BITS;	//18
parameter WINDOW_S = (2*ALPHA+1);

//parameter NUMOFROWS = 1080+ALPHA;	//1080p :^)
//parameter NUMOFCOLS = STRIPEWIDTH+4*ALPHA;	//ex: 180

parameter ROWBITS = $clog2(NUMOFROWS);
parameter COLBITS = $clog2(NUMOFCOLS);

input 					clk;
input 					rst_n;
input 					valid_in;
input signed [I_BITS-1:0] 		Iq;
input signed [I_BITS-1:0] 		Is;

output reg signed [SIGI_BITS-1:0]	SigIi;
output reg 					valid_out;
reg [SIGI_BITS-1:0] SigIi_nxt;

output reg 				 		wen, ren1, ren2;
output reg [COLBITS-1:0] 	 	waddr, raddr1, raddr2;
output reg [II_BITS-1:0] 	 	wdata;
input signed [II_BITS-1:0]	 	rdata1, rdata2;


reg [ROWBITS-1:0] row_count, row_count_nxt;
reg [COLBITS-1:0] col_count, col_count_nxt;

reg signed [II_BITS-1:0] IId_apos, IId_apos_nxt;
reg signed [II_BITS-1:0] IIs;	//IId <= IIs
reg signed [II_BITS-1:0] IId;
reg signed [II_BITS-1:0] IIa_pip, IIa_pip_nxt;


reg 			  valid_in_pip, valid_in_pip_nxt;

always@(*)begin
	row_count_nxt = row_count;
	col_count_nxt = col_count;
	wen = 0; ren1 = 1; ren2 = 1;
	raddr1 = col_count;
	raddr2 = col_count-(WINDOW_S);	
	IIs 		= ((col_count==0)?0:IId)+((row_count==0)?0:rdata1)-((row_count==0||col_count==0)?0:IId_apos)-((row_count<WINDOW_S)?0:Iq)+((row_count>=1080)?0:Is);
	IIa_pip_nxt = ((col_count<WINDOW_S)?0:rdata2);
	SigIi_nxt 	= IId-IIa_pip;	//(pipelined)IIs-IIa
	IId_apos_nxt = rdata1;
	waddr = col_count;
	wdata = IIs;
	valid_in_pip_nxt = (row_count>=ALPHA&&col_count>=(2*ALPHA))?valid_in:0;
	if(valid_in) begin
		row_count_nxt = (col_count==(NUMOFCOLS-1))?((row_count==(NUMOFROWS-1))?0:row_count+1):row_count;
		col_count_nxt = (col_count==(NUMOFCOLS-1))?0:col_count+1;
		raddr1 = col_count_nxt;				// rdata1 : IIs_apos
		raddr2 = col_count_nxt-(WINDOW_S);	// rdata2 : IIa
		wen = 1;
	end
end

always@(posedge clk or negedge rst_n)begin
	if(!rst_n)begin
		row_count <= 0;
		col_count <= 0;
		IId_apos  <= 0;
		IId		  <= 0;
		valid_out <= 0;
		valid_in_pip <= 0;
		SigIi	  <= 0;
		IIa_pip   <= 0;
	end	else begin
		row_count <= row_count_nxt;
		col_count <= col_count_nxt;	
		IId_apos  <= IId_apos_nxt;
		IId		  <= IIs;
		valid_out <= valid_in_pip;
		valid_in_pip <= valid_in_pip_nxt;
		SigIi	  <= SigIi_nxt;
		IIa_pip   <= IIa_pip_nxt;
	end
end


endmodule

//===========================================================================//
//===========================================================================//
module Second_Stage_BK
#(parameter STRIPEWIDTH = 120, parameter ALPHA = 15, parameter A_FRACBITS = 8, parameter B_FRACBITS = 1, parameter A_FRAC_REMAIN = 6)	//S = 2*alpha+1
(clk, rst_n, ak_in, SigIi, Sigpi, valid_in, ak, bk, valid_out);

parameter SIGI_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+8;	//18
parameter AREA_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1));	//10

parameter NUMOFROWS = 1080;	//1080p :^)
parameter NUMOFCOLS = STRIPEWIDTH+2*ALPHA;	//ex: 180

parameter ROWBITS = $clog2(NUMOFROWS);
parameter COLBITS = $clog2(NUMOFCOLS);

parameter NUM_STAGES = 4;

input 				clk;
input 				rst_n;
input 				valid_in;
input [SIGI_BITS-1:0]	SigIi;
input [SIGI_BITS-1:0]	Sigpi;
input signed [4+A_FRACBITS:0]	ak_in;
output valid_out;
output reg signed [4+A_FRACBITS:0]	ak;
output reg signed [12+B_FRACBITS:0]	bk;

wire signed [SIGI_BITS:0]	s_SigIi = {1'b0, SigIi};
wire signed [SIGI_BITS+A_FRAC_REMAIN:0]	s_Sigpi = {1'b0, Sigpi, {A_FRAC_REMAIN{1'b0}}};

reg [ROWBITS-1:0] row_count, row_count_nxt;
reg [COLBITS-1:0] col_count, col_count_nxt;

reg [AREA_BITS-1:0] area_count, area_count_nxt;
wire signed [4+A_FRAC_REMAIN:0] ak_remain = ak_in[4+A_FRACBITS:A_FRACBITS-A_FRAC_REMAIN];
wire signed [4+A_FRAC_REMAIN+SIGI_BITS+1:0] Iiak = s_SigIi*ak_remain;
wire signed [4+A_FRAC_REMAIN+SIGI_BITS+1:0] Numer = s_Sigpi-Iiak;

//wire signed [AREA_BITS+A_FRAC_REMAIN-B_FRACBITS:0] Denom = {1'b0, area_count, {A_FRAC_REMAIN-B_FRACBITS{1'b0}}};
wire signed [AREA_BITS:0] Denom = {1'b0, area_count};

reg signed [4+A_FRAC_REMAIN+SIGI_BITS+1:0] div_a;
reg signed [AREA_BITS:0] div_b;

wire signed [4+A_FRAC_REMAIN+SIGI_BITS+1:0] quotient_inst; //4+3+18+1:0 ->27 bits

wire signed [4+A_FRACBITS:0]	ak_nxt;
reg signed  [12+B_FRACBITS:0]	bk_nxt;

reg div_in_en;

always@(*)begin
	row_count_nxt = row_count;
	col_count_nxt = col_count;
	area_count_nxt = area_count;
	bk_nxt = quotient_inst[12+A_FRAC_REMAIN:A_FRAC_REMAIN-B_FRACBITS];
	if(valid_in) begin
		row_count_nxt = (col_count==(NUMOFCOLS-1))?((row_count==(NUMOFROWS-1))?0:row_count+1):row_count;
		col_count_nxt = (col_count==(NUMOFCOLS-1))?0:col_count+1;
		area_count_nxt = (col_count==(NUMOFCOLS-1))?((row_count<=(ALPHA-1))?area_count+(2*ALPHA+1):((row_count>=(1079-ALPHA))?area_count-(2*ALPHA+1):area_count)):area_count;
		area_count_nxt = (col_count==(NUMOFCOLS-1) && row_count==(NUMOFROWS-1))?(ALPHA+1)*(2*ALPHA+1):area_count_nxt;
	end
end

always@(posedge clk or negedge rst_n)begin
	if(!rst_n)begin
		row_count <= 0;
		col_count <= 0;
		area_count <= (ALPHA+1)*(2*ALPHA+1);
		div_in_en <= 0;
		ak <= 0;
		bk <= 0;
		div_a <= 0;
		div_b <= 0;
	end	else begin
		row_count <= row_count_nxt;
		col_count <= col_count_nxt;	
		area_count <= area_count_nxt;
		div_in_en <= valid_in;
		ak <= ak_nxt;
		bk <= bk_nxt;
		div_a <= Numer;
		div_b <= Denom;
	end
end

//4+6+18+2=30; 10+1= 11
DW_div_pipe #(4+A_FRAC_REMAIN+SIGI_BITS+2, AREA_BITS+1, 1, 1, NUM_STAGES, 1, 1, 0)  U1 
(
.clk(clk),   .rst_n(rst_n),   .en(1'b1),
.a(div_a),   .b(div_b),   .quotient(quotient_inst),
.remainder(),   .divide_by_0() 
);

Shift_Reg #(NUM_STAGES, 1) valid_sr
(.clk(clk), .rst_n(rst_n), .in_data(div_in_en), .out_data(valid_out));

Shift_Reg #(NUM_STAGES, 4+A_FRACBITS+1) ak_sr
(.clk(clk), .rst_n(rst_n), .in_data(ak_in), .out_data(ak_nxt));

endmodule

//==============================================================================//
//===========================================================================//

module Second_Stage_AK
#(parameter STRIPEWIDTH = 120, parameter ALPHA = 15, parameter A_FRACBITS = 8)	//S = 2*alpha+1
(clk, rst_n, SigIi, Sigpi, SigIiIi, SigIipi, valid_in, ak, SigIi_pip, Sigpi_pip, valid_out);

parameter SIGI_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+8;	//18
parameter SIGII_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+16;	//26
parameter AREA_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1));	//10
parameter DENOM_BITS = SIGII_BITS+AREA_BITS+1;	//37
parameter W2EPS = (36'd6005195);//eps*255*255*(2*alpha+1)^4

parameter NUMOFROWS = 1080;	//1080p :^)
parameter NUMOFCOLS = STRIPEWIDTH+2*ALPHA;	//ex: 180

parameter ROWBITS = $clog2(NUMOFROWS);
parameter COLBITS = $clog2(NUMOFCOLS);

parameter NUM_STAGES = 10;

input 				clk;
input 				rst_n;
input 				valid_in;
input [SIGI_BITS-1:0]	SigIi;
input [SIGI_BITS-1:0]	Sigpi;
input [SIGII_BITS-1:0]	SigIiIi;
input [SIGII_BITS-1:0]	SigIipi;
output reg signed [4+A_FRACBITS:0]	ak;
output [SIGI_BITS-1:0]	SigIi_pip;
output [SIGI_BITS-1:0]	Sigpi_pip;
output valid_out;

reg [ROWBITS-1:0] row_count, row_count_nxt;
reg [COLBITS-1:0] col_count, col_count_nxt;

reg [AREA_BITS-1:0] area_count, area_count_nxt;
wire signed [DENOM_BITS-1:0] Ii2Area, IipiArea, IiIi, Iipi;
wire signed [DENOM_BITS-1:0] Denom, Numer;

reg signed [DENOM_BITS+A_FRACBITS-1:0] div_a;
reg signed [DENOM_BITS-1:0] div_b;

wire signed [DENOM_BITS+A_FRACBITS-1:0] quotient_inst;

reg signed [4+A_FRACBITS:0]	ak_nxt;

reg div_in_en;


assign Ii2Area  = area_count*SigIiIi;
assign IipiArea = area_count*SigIipi;
assign IiIi = SigIi*SigIi;
assign Iipi = SigIi*Sigpi;
assign Denom = W2EPS+Ii2Area-IiIi;
assign Numer = IipiArea-Iipi;

always@(*)begin
	row_count_nxt = row_count;
	col_count_nxt = col_count;
	area_count_nxt = area_count;
	ak_nxt = quotient_inst[4+A_FRACBITS:0];
	if(valid_in) begin
		row_count_nxt = (col_count==(NUMOFCOLS-1))?((row_count==(NUMOFROWS-1))?0:row_count+1):row_count;
		col_count_nxt = (col_count==(NUMOFCOLS-1))?0:col_count+1;
		area_count_nxt = (col_count==(NUMOFCOLS-1))?((row_count<=(ALPHA-1))?area_count+(2*ALPHA+1):((row_count>=(1079-ALPHA))?area_count-(2*ALPHA+1):area_count)):area_count;
		area_count_nxt = (col_count==(NUMOFCOLS-1) && row_count==(NUMOFROWS-1))?(ALPHA+1)*(2*ALPHA+1):area_count_nxt;
	end
end

always@(posedge clk or negedge rst_n)begin
	if(!rst_n)begin
		row_count <= 0;
		col_count <= 0;
		area_count <= (ALPHA+1)*(2*ALPHA+1);
		div_in_en <= 0;
		ak <= 0;
		div_a <= 0;
		div_b <= 0;
	end	else begin
		row_count <= row_count_nxt;
		col_count <= col_count_nxt;	
		area_count <= area_count_nxt;
		div_in_en <= valid_in;
		ak <= ak_nxt;
		div_a <= {Numer, {A_FRACBITS{1'b0}}};
		div_b <= Denom;
	end
end

//45/37
DW_div_pipe #(DENOM_BITS+A_FRACBITS, DENOM_BITS, 1, 1, NUM_STAGES, 1, 1, 0)  U1 
(
.clk(clk),   .rst_n(rst_n),   .en(1'b1),
.a(div_a),   .b(div_b),   .quotient(quotient_inst),
.remainder(),   .divide_by_0() 
);

Shift_Reg #(NUM_STAGES, 1) valid_sr
(.clk(clk), .rst_n(rst_n), .in_data(div_in_en), .out_data(valid_out));

Shift_Reg #(NUM_STAGES+1, SIGI_BITS*2) sigIi_sigpi_sr
(.clk(clk), .rst_n(rst_n), .in_data({SigIi, Sigpi}), .out_data({SigIi_pip, Sigpi_pip}));

endmodule

//==============================================================================//
module Shift_Reg#(parameter PIP_NUM = 8, parameter DATA_BITS = 128)
(clk, rst_n, in_data, out_data);

	input 		 clk;
	input 		 rst_n;
	input  [DATA_BITS-1:0]		in_data;
	output [DATA_BITS-1:0]		out_data;

	reg [DATA_BITS-1:0] buf_reg      [PIP_NUM-1:0];
	reg [DATA_BITS-1:0] buf_reg_nxt  [PIP_NUM-1:0];

	integer i;

	assign out_data = buf_reg[PIP_NUM-1];

	always@(posedge clk or negedge rst_n)begin
		if(!rst_n)begin
			for(i=0;i<PIP_NUM;i=i+1)begin
				buf_reg[i] <= 0;
			end
		end
		else begin
			for(i=0;i<PIP_NUM;i=i+1)begin
				buf_reg[i] <= buf_reg_nxt[i];
			end
		end
	end

	always@(*)begin
		buf_reg_nxt[0] = in_data;
		for(i=0;i<PIP_NUM-1;i=i+1)begin
			buf_reg_nxt[i+1] = buf_reg[i];
		end
	end
	
endmodule
//===========================================================================//

module First_Stage
#(parameter STRIPEWIDTH = 120, parameter ALPHA = 15)	//S = 2*alpha+1
(clk, rst_n, Iq, Is, pq, ps, valid_in, SigIi, Sigpi, SigIiIi, SigIipi, valid_out, ready, read_addr_q, read_addr_s);

parameter II_BITS 	= $clog2((2*ALPHA+1)*(STRIPEWIDTH+4*ALPHA))+8;	//21
parameter III_BITS 	= $clog2((2*ALPHA+1)*(STRIPEWIDTH+4*ALPHA))+16;	//29

parameter ADDR_BITS = $clog2(1080*(STRIPEWIDTH+4*ALPHA));	//18
parameter SIGI_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+8;	//18
parameter SIGII_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+16;	//26
parameter WINDOW_S = (2*ALPHA+1);

parameter NUMOFROWS = 1080+ALPHA;	//1080p :^)
parameter NUMOFCOLS = STRIPEWIDTH+4*ALPHA;	//ex: 180

parameter ROWBITS = $clog2(NUMOFROWS);
parameter COLBITS = $clog2(NUMOFCOLS);

input 				clk;
input 				rst_n;
input 				valid_in;
input [7:0] 		Iq;
input [7:0] 		Is;
input [7:0] 		pq;
input [7:0] 		ps;
input				ready; //indicate testbench acknowledge read_addr

//Read Iq/pq
output reg [ADDR_BITS-1:0]  read_addr_q;
//Read Is/ps
output reg [ADDR_BITS-1:0]  read_addr_s;

wire [ADDR_BITS-1:0]  read_addr_s_nxt = (!ready)?read_addr_s:read_addr_s+1;
wire [ADDR_BITS-1:0]  read_addr_q_nxt = (!ready)?read_addr_q:read_addr_s+1-((STRIPEWIDTH+4*ALPHA)*(2*ALPHA+1));

output [SIGI_BITS-1:0]	SigIi;
output [SIGI_BITS-1:0]	Sigpi;
output [SIGII_BITS-1:0]	SigIiIi;
output [SIGII_BITS-1:0]	SigIipi;
output 					valid_out;


wire 				 wen, ren1, ren2;
wire [COLBITS-1:0] 	 waddr, raddr1, raddr2;
wire [II_BITS-1:0] 	 wdata_Ii, wdata_pi;
wire [II_BITS-1:0]	 rdata1_Ii, rdata2_Ii;
wire [II_BITS-1:0]	 rdata1_pi, rdata2_pi;

wire [III_BITS-1:0]	 wdata_IiIi, wdata_Iipi;
wire [III_BITS-1:0]	 rdata1_IiIi, rdata2_IiIi;
wire [III_BITS-1:0]	 rdata1_Iipi, rdata2_Iipi;

reg valid_in_pip;
reg [7:0] Iq_pip, Is_pip, pq_pip, ps_pip;
reg [15:0] IqIq;
reg [15:0] IsIs;
reg [15:0] Iqpq;
reg [15:0] Isps;

//always@(*)begin	
//end


always@(posedge clk or negedge rst_n)begin
	if(!rst_n)begin
		valid_in_pip <= 0;
		Iq_pip <= 0;
		Is_pip <= 0;
		pq_pip <= 0;
		ps_pip <= 0;
		IqIq   <= 0;
		IsIs   <= 0;
		Iqpq   <= 0;
		Isps   <= 0;
		read_addr_s <= 0;
		read_addr_q <= 0;
	end	else begin
		valid_in_pip <= valid_in;
		Iq_pip <= Iq;
		Is_pip <= Is;
		pq_pip <= pq;
		ps_pip <= ps;
		IqIq   <= Iq*Iq;
		IsIs   <= Is*Is;
		Iqpq   <= Iq*pq;
		Isps   <= Is*ps;
		read_addr_s <= read_addr_s_nxt;
		read_addr_q <= read_addr_q_nxt;
	end
end

IIE #(8, ALPHA, II_BITS, SIGI_BITS, NUMOFROWS, NUMOFCOLS) IIE_Ii
(
.clk(clk), .rst_n(rst_n), .Iq(Iq_pip), .Is(Is_pip), .valid_in(valid_in_pip), .SigIi(SigIi), .valid_out(valid_out), 
.wen(wen), .ren1(ren1), .ren2(ren2), .waddr(waddr), .raddr1(raddr1), .raddr2(raddr2), .wdata(wdata_Ii), .rdata1(rdata1_Ii), .rdata2(rdata2_Ii)
);

IIE #(8, ALPHA, II_BITS, SIGI_BITS, NUMOFROWS, NUMOFCOLS) IIE_pi
(
.clk(clk), .rst_n(rst_n), .Iq(pq_pip), .Is(ps_pip), .valid_in(valid_in_pip), .SigIi(Sigpi), .valid_out(), 
.wen(), .ren1(), .ren2(), .waddr(), .raddr1(), .raddr2(), .wdata(wdata_pi), .rdata1(rdata1_pi), .rdata2(rdata2_pi)
);

IIE #(16, ALPHA, III_BITS, SIGII_BITS, NUMOFROWS, NUMOFCOLS) IIE_IiIi
(
.clk(clk), .rst_n(rst_n), .Iq(IqIq), .Is(IsIs), .valid_in(valid_in_pip), .SigIi(SigIiIi), .valid_out(), 
.wen(), .ren1(), .ren2(), .waddr(), .raddr1(), .raddr2(), .wdata(wdata_IiIi), .rdata1(rdata1_IiIi), .rdata2(rdata2_IiIi)
);

IIE #(16, ALPHA, III_BITS, SIGII_BITS, NUMOFROWS, NUMOFCOLS) IIE_Iipi
(
.clk(clk), .rst_n(rst_n), .Iq(Iqpq), .Is(Isps), .valid_in(valid_in_pip), .SigIi(SigIipi), .valid_out(), 
.wen(), .ren1(), .ren2(), .waddr(), .raddr1(), .raddr2(), .wdata(wdata_Iipi), .rdata1(rdata1_Iipi), .rdata2(rdata2_Iipi)
);

//21+29+21+29=100, 180 words
memory_wrapper  #(II_BITS+II_BITS+III_BITS+III_BITS, NUMOFCOLS) rf_mem
(
	.clk(clk), .rst_n(rst_n), 
	.wen(wen), .waddr(waddr), .wdata({wdata_Ii, wdata_pi, wdata_IiIi, wdata_Iipi}), 
	.ren1(ren1), .raddr1(raddr1), .rdata1({rdata1_Ii, rdata1_pi, rdata1_IiIi, rdata1_Iipi}), 
	.ren2(ren2), .raddr2(raddr2), .rdata2({rdata2_Ii, rdata2_pi, rdata2_IiIi, rdata2_Iipi})
);

endmodule

//===========================================================================//
module IIE
#(parameter I_BITS = 8, parameter ALPHA = 15, parameter II_BITS = 21, parameter SIGI_BITS = 18, parameter NUMOFROWS =1095, parameter NUMOFCOLS = 180)	//S = 2*alpha+1
(clk, rst_n, Iq, Is, valid_in, SigIi, valid_out, wen, ren1, ren2, waddr, raddr1, raddr2, wdata, rdata1, rdata2);

//parameter II_BITS 	= $clog2((2*ALPHA+1)*(STRIPEWIDTH+4*ALPHA))+I_BITS;	//21
//parameter SIGI_BITS	= $clog2((2*ALPHA+1)*(2*ALPHA+1))+I_BITS;	//18
parameter WINDOW_S = (2*ALPHA+1);

//parameter NUMOFROWS = 1080+ALPHA;	//1080p :^)
//parameter NUMOFCOLS = STRIPEWIDTH+4*ALPHA;	//ex: 180

parameter ROWBITS = $clog2(NUMOFROWS);
parameter COLBITS = $clog2(NUMOFCOLS);

input 					clk;
input 					rst_n;
input 					valid_in;
input [I_BITS-1:0] 		Iq;
input [I_BITS-1:0] 		Is;

output reg [SIGI_BITS-1:0]	SigIi;
output reg 					valid_out;
reg [SIGI_BITS-1:0] SigIi_nxt;

output reg 				 wen, ren1, ren2;
output reg [COLBITS-1:0] 	 waddr, raddr1, raddr2;
output reg [II_BITS-1:0] 	 wdata;
input [II_BITS-1:0]	 rdata1, rdata2;


reg [ROWBITS-1:0] row_count, row_count_nxt;
reg [COLBITS-1:0] col_count, col_count_nxt;

reg [II_BITS-1:0] IId_apos, IId_apos_nxt;
reg [II_BITS-1:0] IIs;	//IId <= IIs
reg [II_BITS-1:0] IId;
reg [II_BITS-1:0] IIa_pip, IIa_pip_nxt;


reg 			  valid_in_pip, valid_in_pip_nxt;

always@(*)begin
	row_count_nxt = row_count;
	col_count_nxt = col_count;
	wen = 0; ren1 = 1; ren2 = 1;
	raddr1 = col_count;
	raddr2 = col_count-(WINDOW_S);	
	IIs 		= ((col_count==0)?0:IId)+((row_count==0)?0:rdata1)-((row_count==0||col_count==0)?0:IId_apos)-((row_count<WINDOW_S)?0:Iq)+((row_count>=1080)?0:Is);
	IIa_pip_nxt = ((col_count<WINDOW_S)?0:rdata2);
	SigIi_nxt 	= IId-IIa_pip;	//(pipelined)IIs-IIa
	IId_apos_nxt = rdata1;
	waddr = col_count;
	wdata = IIs;
	valid_in_pip_nxt = (row_count>=ALPHA&&col_count>=(2*ALPHA))?valid_in:0;
	if(valid_in) begin
		row_count_nxt = (col_count==(NUMOFCOLS-1))?((row_count==(NUMOFROWS-1))?0:row_count+1):row_count;
		col_count_nxt = (col_count==(NUMOFCOLS-1))?0:col_count+1;
		raddr1 = col_count_nxt;				// rdata1 : IIs_apos
		raddr2 = col_count_nxt-(WINDOW_S);	// rdata2 : IIa
		wen = 1;
	end
end

always@(posedge clk or negedge rst_n)begin
	if(!rst_n)begin
		row_count <= 0;
		col_count <= 0;
		IId_apos  <= 0;
		IId		  <= 0;
		valid_out <= 0;
		valid_in_pip <= 0;
		SigIi	  <= 0;
		IIa_pip   <= 0;
	end	else begin
		row_count <= row_count_nxt;
		col_count <= col_count_nxt;	
		IId_apos  <= IId_apos_nxt;
		IId		  <= IIs;
		valid_out <= valid_in_pip;
		valid_in_pip <= valid_in_pip_nxt;
		SigIi	  <= SigIi_nxt;
		IIa_pip   <= IIa_pip_nxt;
	end
end


endmodule

//===========================================================================//
//raddr1 and raddr2 have to be one odd, one even
module memory_wrapper #(parameter DATA_WIDTH=100, parameter WORD_WIDTH=180)
(clk, rst_n, wen, waddr, wdata, ren1, raddr1, rdata1, ren2, raddr2, rdata2);

parameter ADDR_WIDTH = $clog2(WORD_WIDTH)-1;

input wen, ren1, ren2, clk, rst_n;
input [ADDR_WIDTH:0] waddr, raddr1, raddr2;
input [(DATA_WIDTH-1):0] wdata;
output [(DATA_WIDTH-1):0] rdata1, rdata2;
	
wire ren_even = (ren1&&(!raddr1[0])) || (ren2&&(!raddr2[0]));
wire ren_odd  = (ren1&&(raddr1[0]))  || (ren2&&(raddr2[0]));
wire [ADDR_WIDTH-1:0] raddr_even = (ren1&&(!raddr1[0]))?(raddr1[(ADDR_WIDTH):1]):(raddr2[(ADDR_WIDTH):1]);
wire [ADDR_WIDTH-1:0] raddr_odd  = (ren1&&(raddr1[0]))? (raddr1[(ADDR_WIDTH):1]):(raddr2[(ADDR_WIDTH):1]);

wire [(DATA_WIDTH-1):0] rdata_even, rdata_odd;

//pipelined memory address to match output data
reg raddr1_is_odd;

assign rdata1 = (raddr1_is_odd)?rdata_odd:rdata_even;
assign rdata2 = (raddr1_is_odd)?rdata_even:rdata_odd;

wire wen_even = (wen&&(!waddr[0]));
wire wen_odd  = (wen&&(waddr[0]));


always@(posedge clk or negedge rst_n)begin
	if(!rst_n)begin
		raddr1_is_odd <= 0;
	end
	else begin
		raddr1_is_odd <= raddr1[0];
	end
end


data_mem  #(DATA_WIDTH, WORD_WIDTH>>1) register_file_even	//even addr (0, 2, 4, 6)
(.clk(clk), .rst_n(rst_n), .wen(wen_even), .waddr(waddr[(ADDR_WIDTH):1]), .wdata(wdata), .ren(ren_even), .raddr(raddr_even), .rdata(rdata_even));

data_mem  #(DATA_WIDTH, WORD_WIDTH>>1) register_file_odd //odd addr (1, 3, 5, 7)
(.clk(clk), .rst_n(rst_n), .wen(wen_odd), .waddr(waddr[(ADDR_WIDTH):1]), .wdata(wdata), .ren(ren_odd), .raddr(raddr_odd), .rdata(rdata_odd));

endmodule

//===========================================================================//


module data_mem
#(parameter DATA_WIDTH=100, parameter WORD_WIDTH=180)
(clk, rst_n, wen, waddr, wdata, ren, raddr, rdata);

	parameter ADDR_WIDTH = $clog2(WORD_WIDTH);	//9

	input wen, ren, clk, rst_n;
	input [(ADDR_WIDTH-1):0] waddr, raddr;
	input [(DATA_WIDTH-1):0] wdata;
	output reg [(DATA_WIDTH-1):0] rdata;
	
	integer i;
	
	// Declare the RAM variable
	reg [DATA_WIDTH-1:0] ram[WORD_WIDTH-1:0];
	// Port A 
	always @ (posedge clk or negedge rst_n)
	begin
		if(!rst_n) begin
			for(i=0; i< WORD_WIDTH; i= i+1)begin
				ram[i] <= 0;
			end
		end else begin
			if (wen) begin
				ram[waddr] <= wdata;
			end 
			if (ren) begin
				rdata <= ram[raddr];
			end 
		end
	end 
endmodule 