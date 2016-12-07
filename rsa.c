#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include <stdio.h>
#include <time.h>

using namespace std;
int primetest(mpz_t p1);
void euclides(mpz_t result,mpz_t a, mpz_t b);
void chinese(mpz_t res,mpz_t p1,mpz_t p2,mpz_t uz,mpz_t pow);
void generate_random(mpz_t random1,mpz_t random2,mpz_t random3);
void fastpow(mpz_t res,mpz_t base,mpz_t pow,mpz_t mod);

int main() {
    mpz_t r1,r2,r3,pmin1,ppmin1,fi,one,sk,omsg,pmsg,n,pk,res,two;
    mpz_init(fi);
    mpz_init(res);
    mpz_init(pmin1);
    mpz_init(ppmin1);
    mpz_init(one);
    mpz_init(omsg);
    mpz_init(pmsg);
    mpz_init(n);
    mpz_init(pk);
    mpz_init(sk);
    mpz_init_set_str(one,"1",10);
    mpz_init_set_str(two,"2",10);
    mpz_init(r1);
    mpz_init(r2);
    mpz_init(r3);
    generate_random(r1,r2,r3);

    do {
        if(mpz_even_p(r1)) {
            mpz_sub(r1,r1,one);
        }
        if(mpz_even_p(r2)) {
            mpz_sub(r2,r2,one);
        }
        if(mpz_even_p(r3)) {
            mpz_sub(r3,r3,one);
        }
        if(!(primetest(r1)))
            mpz_add(r1,r1,two);
        if(!(primetest(r2)))
            mpz_add(r2,r2,two);
        if(!(primetest(r3)))
            mpz_add(r3,r3,two);
    } while(!(primetest(r1)&&(primetest(r2)&&(primetest(r3)))));

    mpz_mul(n,r1,r2);
    mpz_sub(pmin1,r1,one);
    mpz_sub(ppmin1,r2,one);
    mpz_mul(fi,pmin1,ppmin1);
    mpz_set(sk,r3);
    euclides(pk,fi,r3);

    gmp_printf("Public key %Zd\n",pk);
    do {
        gmp_printf("Enter the message: ");
        gmp_scanf("%Zd",pmsg);
    }
    while((mpz_cmp(n,pmsg))<0);
    fastpow(omsg,pmsg,pk,n);
    gmp_printf("Public message: %Zd\n", omsg);
    fastpow(pmsg,omsg,sk,n);
    gmp_printf("Private message: %Zd\n",pmsg);
    chinese(res,r1,r2,omsg,sk);
    gmp_printf("Chinese message : %Zd\n",res);
    return 0;
}
// MILLER RABIN PRIMETEST
int primetest(mpz_t p1) {
    mpz_t two,one,zero,evennum,oddnum,pow,p1_min1;
    mpz_init(evennum);
    mpz_init(p1_min1);
    mpz_init(oddnum);
    mpz_init(pow);
    mpz_init_set_str(two,"2",10);
    mpz_init_set_str(one,"1",10);
    mpz_init_set_str(zero,"0",10);
    mpz_sub_ui(p1_min1,p1,1);
    mpz_set(evennum,p1_min1);

    while(mpz_even_p(evennum)) {
        mpz_fdiv_q_2exp(evennum,evennum,1);
        mpz_add(pow,pow,one);
        gmp_printf("%Zd\n",evennum);
    }
    mpz_set(oddnum,evennum);
    fastpow(evennum,two,oddnum,p1);
    if(mpz_cmp_ui(evennum,1)==0) {
        return 1;
    }
    for(; mpz_cmp(pow,zero)>0; mpz_sub(pow,pow,one)) {
        if(mpz_cmp(evennum,p1_min1)==0) {
            return 1;
        }
        fastpow(evennum,evennum,two,p1);
    }
    if(mpz_cmp(evennum,p1_min1)==0) {
        return 1;
    }
    return 0;
}

void euclides(mpz_t result,mpz_t a, mpz_t b) {
    mpz_t firsta,q,x,lastx,y,lasty,temp1,temp2,temp3,one,helper,null,a1,b1;
    mpz_init(firsta);
    mpz_init(q);
    mpz_init(x);
    mpz_init(lastx);
    mpz_init(y);
    mpz_init(lasty);
    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init(temp3);
    mpz_init(one);
    mpz_init(helper);
    mpz_init(null);
    mpz_init(a1);
    mpz_init(b1);
    mpz_set(a1,a);
    mpz_set(b1,b);
    mpz_set(firsta,a);
    mpz_init_set_str(x,"0",10);
    mpz_init_set_str(y,"1",10);
    mpz_init_set_str(lastx,"1",10);
    mpz_init_set_str(lasty,"0",10);
    mpz_init_set_str(one,"1",10);
    mpz_init_set_str(null,"0",10);

    while (mpz_cmp(b1,one)!=0) {
        mpz_fdiv_q(q,a1,b1);
        mpz_mod(temp1,a1,b1);
        mpz_set(a1,b1);
        mpz_set(b1,temp1);

        mpz_set(temp2,x);
        mpz_mul(helper,q,x);
        mpz_sub(x,lastx,helper);
        mpz_set(lastx,temp2);

        mpz_set(temp3,y);
        mpz_mul(helper,q,y);

        mpz_sub(y,lasty,helper);
        mpz_set(lasty,temp3);
    }
    //IF Y IS POSITIVE, Y=Y, IF NEGATIVE, THEN A+Y
    if(mpz_cmp(y,null)) {
        mpz_add(helper,firsta,y);
        mpz_set(result,helper);
    } else {
        mpz_set(result,y);
    }
}

void chinese(mpz_t res,mpz_t p1,mpz_t p2,mpz_t uz,mpz_t pow) {
    mpz_t y1,y2,min1,a1,a2,x1,x2,x,temp1,temp2,pp1;
    mpz_init(a1);
    mpz_init(x1);
    mpz_init(x2);
    mpz_init(x);
    mpz_init(a2);
    mpz_init(y1);
    mpz_init(y2);
    mpz_init(pp1);
    mpz_init(temp1);
    mpz_init(temp2);
    mpz_init(min1);
    mpz_init_set_ui(min1,-1);
    fastpow(y1,p1,min1,p2);
    fastpow(y2,p2,min1,p1);
    mpz_mul(pp1,p1,p2);
    fastpow(a1,uz,pow,p2);
    fastpow(a2,uz,pow,p1);
    mpz_mul(temp1,a1,p1);
    mpz_mul(temp2,temp1,y1);
    mpz_set(x1,temp2);
    mpz_mul(temp1,a2,p2);
    mpz_mul(temp2,temp1,y2);
    mpz_add(x,x1,temp2);
    mpz_mod(res,x,pp1);
}

void fastpow(mpz_t res,mpz_t base,mpz_t pow,mpz_t mod) {
    mpz_t b,zero,temp,temp2,temp3,gotpow,one,two,base,modulo,min1;
    mpz_init_set_ui(b,1);
    mpz_init_set_ui(one,1);
    mpz_init_set_ui(two,2);
    mpz_init_set_ui(zero,0);
    mpz_init(min1);
    mpz_sub(min1,zero,one);
    mpz_init_set(gotpow,pow);
    mpz_init_set(base,base);
    mpz_init_set(modulo,mod);
    mpz_init(temp);
    mpz_init(temp2);
    mpz_init(temp3);
    mpz_set(temp,base);

    if((mpz_cmp(gotpow,zero)==0)&&(mpz_cmp(base,zero)==0)) {
        printf("Nem meghatarozhato");
    }
    if(mpz_cmp(gotpow,zero)==0) {
        mpz_set(res,b);
    }

    mpz_mod(temp3,gotpow,two);
    if(mpz_cmp(temp3,one)==0)
        mpz_set(b,base);
    mpz_tdiv_q_2exp(gotpow,gotpow,1);
    while(mpz_cmp(gotpow,zero)!=0) {
        mpz_mul(temp,temp,temp);
        mpz_mod(temp,temp,modulo);
        mpz_mod(temp2,gotpow,two);
        if(mpz_cmp(temp2,one)==0) {
            mpz_mul(b,temp,b);
            mpz_mod(b,b,modulo);
        }
        mpz_tdiv_q_2exp(gotpow,gotpow,1);
    }
    mpz_set(res,b);

    if(mpz_cmp_ui(pow,-1)==0) {
        mpz_invert(res,base,mod);
        mpz_mod(res,res,mod);
    }
}

void generate_random(mpz_t random1,mpz_t random2,mpz_t random3) {
    srand(time(NULL));
    mpz_t rand_Num;
    mpz_t rand_Num2;
    mpz_t rand_Num3;
    unsigned long int seed;
    gmp_randstate_t r_state;

    seed = rand();
    gmp_randinit_mt (r_state);
    gmp_randseed_ui(r_state, seed);

    mpz_init(rand_Num);
    mpz_init(rand_Num2);
    mpz_init(rand_Num3);
    mpz_urandomb(rand_Num,r_state,512);
    mpz_urandomb(rand_Num2,r_state,512);
    mpz_urandomb(rand_Num3,r_state,512);

    mpz_set(random1,rand_Num);
    mpz_set(random2,rand_Num2);
    mpz_set(random3,rand_Num3);

    gmp_randclear(r_state);
    mpz_clear(rand_Num);
}
