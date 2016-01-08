#include <Integer.h>
#include <time.h>

#define yes	1
#define no	0

#define DISC(a,b) (Integer) ((4*a*a*a)+(27*b*b)) // biet so cua duong cong 4a^3+27b^2

#define K_incr 1	// buoc nhay cua K
#define a_incr 30	// buoc nhay cua a

#define max_tries_a  40	// So lan thay doi a truoc khi thay doi K

// Khoi tao (x, y) = (1 ,1)  
#define START_X  1	
#define START_Y  1	

//Tinh K tu n ; K=5(log n)^2
Integer optimal_K(Integer n) {
	Integer K;
	int t;
	K = logn(n, 10);
	K = 5*K*K;
	if (K<2) K=2;
	return K;
}

// k= lcm[1,2,3,4,..,K]
Integer compute_k(Integer K) {
	return LCM(K);
}

Integer Factor(Integer n, Integer K) {
    Integer a, b, k, g, x, y, tempx, tempy, sum_x, sum_y, grad;
    int i, bits, newcurve, a_tries, count=0; 

    a = 1;	// khoi tao he so a cua duong cong y^2= x^3 + ax + b
    k = compute_k(K);
    bits = log2(k); 
    a_tries = 0;
    
    for (;;) {
	x = START_X; 
	y = START_Y;	    
        b = (y*y)-(x*x*x)-(a*x);// Tinh he so b cua duong cong

		g = gcd( DISC(a,b), n);	//   g=  gcd(4a^3+27b^2 ,n )
		while (g != 1) {			
		    if (g == n) {		 // neu g=n tang a va tinh lai b , g	
				a += a_incr;
				b = y*y - x*x*x - a*x;
				g = gcd( DISC(a,b), n);
				continue;
		    }
	    	return g;		
		}			             

	// Them diem vao duong cong
	// Tinh  2(x,y), 4(x,y), 8(x,y), ...  va cong vao (sum_x, sum_y)
	
		sum_x = 0; 
		sum_y = 0; 
		newcurve = no;  // newcurve = yes  Can thay doi duong cong voi k moi

		if (a_tries > max_tries_a) {
			newcurve = yes; // Thay doi duong cong khi so lan thay doi a vuot gioi han
		}
		
		else 
			for (i=1; i<=bits; i++) {		// tinh k*(x, y) 
		    	tempx = x; 
				tempy = 2*y;
		   		g = gcd(tempy, n);
		   		if (1<g && g<n) return g;
		   		
		    	if (g == n) {
					newcurve = yes;  // g=n thay doi duong cong voi k moi
					break;
				}
				
				// g=1
		   		grad = InvModN(tempy, n) * ( 3*(x*x) + a); // grad = (tempy)^-1(mod n) * (3x^2+a)
		   		x = (grad*grad - 2*x) % n;             // x = (grad^2 - 2x) mod n
		   		y = (grad*(tempx-x) - y) % n;          // y= (grad(tempx-x) -y) mod n
	
		    	if (!testbit(k, i)) continue;
		    	
		   		if (!sum_y) {
					sum_x = x;    
					sum_y = y; 
					continue;
				}

				//  Cong {x,y} vao {sumx,sumy}
				
			    tempx = (x - sum_x) % n;
			    tempy = (y - sum_y) % n;
			    
			    g = gcd(tempx, n);
			    
			    if (1<g && g<n) return g;
			    
			    if (g == n) {
					newcurve = yes;  // Thay doi duong cong voi k moi
					break;
				}
				
			    grad = InvModN(tempx ,n) * tempy;		// grad =(tempx)^-1(mod n) *temp y 
			    tempx = ( (grad*grad - sum_x) - x)% n;  //tempx = (grad^2 - sumx -x) mod n
			    tempy = ( grad*(sum_x - tempx) - sum_y) % n; // tempy= (grad(sumx-tempx) -sumy) mod n
			    sum_x = tempx; 
				sum_y = tempy;
			}
	
			if (newcurve) {		//Thay doi duong cong voi k moi
				a = 1;				// gan lai a=1, a_tries =0
				a_tries = 0;			
				do {K += K_incr;}   
				while (k == compute_k(K)); // tinh lai k = lcm[1,2,...,K] Khi chua co gia tri k moi tang K
				k = compute_k(K);		
				bits = log2(k);
			}
			
			else {			
				a += a_incr;	 //Tiep tuc vong lap voi a moi
				a_tries++;	
			}
    }
}


int main(int argc, char *argv[]) {
    Integer n, K, p;
    clock_t start, end;
    int t;
    if (argc < 2) {
		cerr << "Usage: " << argv[0] << " NumberToFactor" << endl;
		return -1;
    }

    n = argv[1];
    if (argc == 3) K = argv[2];
    else K = optimal_K(n);		// default neu k khong duoc chon

   
    start = clock();
	p = Factor(n, K);
	end = clock();
	cout << n << " = " << flush;
	cout << p << " x " << flush;
	cout << n/p << "\n" << flush;
    cout <<"Time = "<< (end-start)/CLOCKS_PER_SEC << " s\n" ;
  
    return 0;
}
