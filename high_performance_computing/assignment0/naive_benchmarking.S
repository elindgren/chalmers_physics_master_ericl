	.file	"naive_benchmarking.c"
	.text
	.section	.rodata.str1.1,"aMS",@progbits,1
.LC0:
	.string	"t=%f\n"
	.section	.rodata.str1.8,"aMS",@progbits,1
	.align 8
.LC3:
	.string	"Elapsed time: %fs - Time per iteration: %fs \n"
	.section	.text.startup,"ax",@progbits
	.p2align 4
	.globl	main
	.type	main, @function
main:
.LFB11:
	.cfi_startproc
	pushq	%rbx
	.cfi_def_cfa_offset 16
	.cfi_offset 3, -16
	subq	$16, %rsp
	.cfi_def_cfa_offset 32
	call	clock
	movq	%rax, %rbx
	call	clock
	pxor	%xmm0, %xmm0
	movl	$.LC0, %edi
	subq	%rbx, %rax
	cvtsi2sdq	%rax, %xmm0
	movl	$1, %eax
	movsd	%xmm0, 8(%rsp)
	call	printf
	movl	$.LC3, %edi
	movl	$2, %eax
	movsd	8(%rsp), %xmm0
	divsd	.LC1(%rip), %xmm0
	movapd	%xmm0, %xmm1
	divsd	.LC2(%rip), %xmm1
	call	printf
	addq	$16, %rsp
	.cfi_def_cfa_offset 16
	xorl	%eax, %eax
	popq	%rbx
	.cfi_def_cfa_offset 8
	ret
	.cfi_endproc
.LFE11:
	.size	main, .-main
	.section	.rodata.cst8,"aM",@progbits,8
	.align 8
.LC1:
	.long	0
	.long	1093567616
	.align 8
.LC2:
	.long	0
	.long	1076101120
	.ident	"GCC: (GNU) 10.1.1 20200507 (Red Hat 10.1.1-1)"
	.section	.note.GNU-stack,"",@progbits
