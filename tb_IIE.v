//behavior tb
`timescale 1ns/10ps
`define CYCLE	10
`define End_CYCLE  10000000            // Modify cycle times once your design need more cycle times!
`define PATI        "./I.dat" 
`define PATP        "./p.dat" 
`define PATQ        "./q.dat"
`define PATIMAGE    "./outputImage_stripe_2"
`define SDFFILE     "./guide_image_syn.sdf"


module tb_IIE;
	parameter Tolerance  = 4'd1;

	parameter STRIPEWIDTH  = 120;
	parameter ALPHA = 15;
	parameter NUMOFCOLS = STRIPEWIDTH+4*ALPHA;	//ex: 180
	parameter NUMOFROWS = 1080+ALPHA;	//1080p :^)
	parameter NUM_PAT_In  = 1080*NUMOFCOLS;
	parameter NUM_PAT_Out = 1080*(STRIPEWIDTH);//+2*ALPHA
	
	parameter A_FRACBITS = 8;
	parameter B_FRACBITS = 1;
	
	//First Stage
	parameter SIGI_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+8;	//18
	parameter SIGII_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+16;	//26
	//Second Stage
	parameter ADDR_BITS = $clog2(1080*(STRIPEWIDTH+4*ALPHA));	//18
	//Third Stage
	parameter ADDR_BITS_THIRD 	= $clog2((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA));	//13
	parameter SIGIA_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+(4+A_FRACBITS+1);	//10+13=23
	parameter SIGIB_BITS = $clog2((2*ALPHA+1)*(2*ALPHA+1))+(12+B_FRACBITS+1);	//10+14=24
	
	
    reg  		 	clk;
	reg  		 	reset;
    reg  [7:0]  	Iq;
    reg  [7:0]  	Is;
    reg  [7:0]  	pq;
    reg  [7:0]  	ps;
    reg  		 	valid_in;
	wire signed [4+A_FRACBITS:0]	ak;
	wire signed [12+B_FRACBITS:0]	bk;
	wire					valid_out, valid_inter_sigab;
	reg 					ready;
	//Read Iq/pq, Is/ps
	wire [ADDR_BITS-1:0]  read_addr_q;
	wire [ADDR_BITS-1:0]  read_addr_s;
	wire [ADDR_BITS-1:0]  read_addr_Iin;
	
	//Third Stage
	wire signed [4+A_FRACBITS:0]	aq;
	wire signed [12+B_FRACBITS:0]	bq;
	wire [ADDR_BITS_THIRD-1:0]		read_addr_q_ab;
	wire signed [SIGIA_BITS-1:0]	Sigak;
	wire signed [SIGIB_BITS-1:0]	Sigbk;
	
	reg  [7:0] 	I_in;
	wire [7:0]	qi;
	
	
	
	parameter ROWBITS = $clog2(NUMOFROWS);
	parameter COLBITS = $clog2(NUMOFCOLS);

	reg [ROWBITS-1:0] row_count;
	reg [COLBITS-1:0] col_count;
	
	reg  [7:0]  			golden_I_in  [NUM_PAT_In-1:0];
	reg  [7:0]  			golden_P_in  [NUM_PAT_In-1:0];
	//reg  signed [SIGIA_BITS-1:0]  	golden_ak_out [NUM_PAT_Out-1:0];
	//reg  signed [SIGIB_BITS-1:0]  	golden_bk_out [NUM_PAT_Out-1:0];
	reg  [7:0]  	golden_q_out [NUM_PAT_Out-1:0];
	
	//Act as temporary Ak/ Bk off-chip memory
	data_mem_imre #(5+A_FRACBITS+13+B_FRACBITS, ((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA))) akbk_mem
    (
	.clk(clk) ,.rst_n(reset), 
	.wen(valid_inter_sigab), 	.waddr(read_addr_q_ab), .wdata({ak,bk}), 
	.ren(1'b1), 				.raddr(read_addr_q_ab), .rdata({aq,bq})
	);
	//wire [5+A_FRACBITS+13+B_FRACBITS-1:0] akbk_temp_mem [((2*ALPHA+1)*(STRIPEWIDTH+2*ALPHA))-1:0];
	
	reg  [31:0]		in_index, out_index;
	
	initial begin
		$display("--------------------------- [ Simulation Starts !! ] ---------------------------");	
		$readmemb(`PATI, golden_I_in);
		$readmemb(`PATP, golden_P_in);
		$readmemb(`PATQ, golden_q_out);
	end
        
	
    `ifdef SDF
	Guided_Image guide_image_u
	`else
	Guided_Image #(STRIPEWIDTH, ALPHA) guide_image_u
	`endif
	(.clk(clk), .rst_n(reset), .valid_in(valid_in), .Iq(Iq), .Is(Is), .pq(pq), .ps(ps), 
	.ak(ak), .bk(bk), .valid_inter_sigab(valid_inter_sigab), .ready(ready), .read_addr_q(read_addr_q), .read_addr_s(read_addr_s),
	.aq(aq), .bq(bq), .Sigak(Sigak), .Sigbk(Sigbk), .read_addr_q_ab(read_addr_q_ab),
	.I_in(I_in), .qi(qi), .read_addr_Iin(read_addr_Iin), .valid_out(valid_out)
	);

	
	`ifdef SDF
	   initial $sdf_annotate(`SDFFILE, guide_image_u );
	`endif
	
    // waveform dump
    initial begin
        $fsdbDumpfile( "guide_image_u.fsdb" );
        $fsdbDumpvars(0, guide_image_u, "+mda");
    end
	
    // clock
    initial begin
        clk = 1'b0;
        forever #(`CYCLE/2) clk = ~clk;
    end
	
	wire stop_all = (out_index==NUM_PAT_Out-1);
	integer err_cnt, i, k, iters, last_valid;
	
	integer im_w;
	
    initial begin
		im_w = $fopen(`PATIMAGE, "w");
		
        reset = 1'b1; valid_in = 1'b0; in_index = 0; out_index = 0; err_cnt = 0; iters = 0; last_valid=1; ready = 1'b0;
        #(`CYCLE) reset = 1'b0;
		#(3*`CYCLE) reset = 1'b1;
		#(3*`CYCLE)
		while(!stop_all) begin
			@(negedge clk);
			iters = iters+1;
			k = $random(iters)%20;
			if(in_index<=((NUMOFROWS*NUMOFCOLS)-1)) begin	//send LDUX
				ready = 1'b1;
				valid_in = 1'b1;
				Iq = golden_I_in[read_addr_q];
				Is = golden_I_in[read_addr_s];
				pq = golden_P_in[read_addr_q];
				ps = golden_P_in[read_addr_s];
				in_index = in_index+1;
			end
			I_in = golden_I_in[read_addr_Iin];
			if(valid_out) begin
				$fwrite(im_w, "%d\n", qi);
				if (!(golden_q_out[out_index]===0&&qi<=Tolerance || golden_q_out[out_index]!==0&&(golden_q_out[out_index]-qi<=Tolerance || qi-golden_q_out[out_index]<=Tolerance))) begin
					$display("Q Pattern %d failed !. Expected candidate = %d, but the Response candidate = %d !! ", out_index, golden_q_out[out_index], qi);
					err_cnt = err_cnt + 1;
				end else begin
					$display("Q Pattern %d is passed !. Expected candidate = %d, Response candidate = %d !! ", out_index, golden_q_out[out_index], qi);
				end
				out_index = out_index+1;
			end
        end
		
		$fclose(im_w);
		
		#(`CYCLE*2); 
		$display("--------------------------- Simulation Stops !!---------------------------");
		if (err_cnt) begin 
			$display("============================================================================");
			$display("\n (T_T) ERROR found!! There are %d errors in total.\n", err_cnt);
			$display("============================================================================");
		end
		 else begin 
			$display("============================================================================");
			$display("\n");
			$display("        ****************************              ");
			$display("        **                        **        /|__/|");
			$display("        **  Congratulations !!    **      / O,O  |");
			$display("        **                        **    /_____   |");
			$display("        **  Simulation Complete!! **   /^ ^ ^ \\  |");
			$display("        **                        **  |^ ^ ^ ^ |w|");
			$display("        *************** ************   \\m___m__|_|");
			$display("\n");
			$display("============================================================================");
			$finish;
		end
		$finish;
    end
	

	always@(err_cnt) begin
		if (err_cnt >= 1) begin
			$display("============================================================================");
			$display("\n (>_<) ERROR!! There are more than 10 errors during the simulation! Please check your code @@ \n");
			$display("============================================================================");
			$finish;
		end
	end
	
	initial begin 
		#`End_CYCLE;
		$display("================================================================================================================");
		$display("(/`n`)/ ~#  There is something wrong with your code!!"); 
		$display("Time out!! The simulation didn't finish after %d cycles!!, Please check it!!!", `End_CYCLE); 
		$display("================================================================================================================");
		$finish;
	end
endmodule

//immediate read/ delay-1 write
module data_mem_imre
#(parameter DATA_WIDTH=100, parameter WORD_WIDTH=180)
(clk, rst_n, wen, waddr, wdata, ren, raddr, rdata);

	parameter ADDR_WIDTH = $clog2(WORD_WIDTH);	//9

	input wen, ren, clk, rst_n;
	input [(ADDR_WIDTH-1):0] waddr, raddr;
	input [(DATA_WIDTH-1):0] wdata;
	output [(DATA_WIDTH-1):0] rdata;
	integer i;
	
	// Declare the RAM variable
	reg [DATA_WIDTH-1:0] ram[WORD_WIDTH-1:0];
	
	assign rdata = (ren)?ram[raddr]:0;
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
		end
	end 
endmodule