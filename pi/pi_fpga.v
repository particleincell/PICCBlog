module pi_fpga;
  reg clk;
  integer i,j;
  parameter N = 800;
  real x,y;
  real fi, fj;
  integer N_in;

  initial begin
    clk <= 0;
    i <= -1; j<= 0;
    N_in = 0;
    // $monitor ("T=%0t clk=%0d, i=%0d, j=%0d, x=%0f", $time, clk, i, j,x);
  end
 
  // create a clock
  always #1 clk = ~clk;
 
  // code to run on each clock "tick"
  always @ (posedge clk) begin
    i = i + 1;
    if (i>=N) begin
      i<=0;     
      j=j+1;
      if (j>=N) begin
        fi = N_in;
        $display("Using %0dx%0d mesh (%0d points), pi = %0f",N,N,N*N,4*fi/(N*N));
        $finish;
      end
    end
    fi = i;
    fj = j;
    x = fi/N;
    y = fj/N;
    if ((x*x+y*y)<=1) N_in <= N_in+1;
  end
endmodule
