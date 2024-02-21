#!/bin/bash
#define directory
Ref=refdata-gex-mm10-2020-A
fas_dir=rawdata
image_dir=HE
slide_dir=slides
res_dir=result

echo 'Start 24AD'
spaceranger count --id=AD335 \
--transcriptome=${Ref} \
--fastqs=${fas_dir}/AD335 \
--image=${image_dir}/AD335.tif \
--sample=AD335 --slide=V11F09-071 \
--slidefile=${slide_dir}/V11F09-071.gpr \
--area=C1 --localcores=2
echo 'AD335 complete'

echo 'Start 24SL'
spaceranger count --id=SL298 \
--transcriptome=${Ref} \
--fastqs=${fas_dir}/SL298 \
--image=${image_dir}/SL298.tif \
--sample=SL298 --slide=V11F09-071 \
--slidefile=${slide_dir}/V11F09-071.gpr \
--area=D1 --localcores=2
echo 'SL298 complete'

echo 'Start 24WT'
spaceranger count --id=WT281 \
--transcriptome=${Ref} \
--fastqs=${fas_dir}/WT281 \
--image=${image_dir}/WT281.tif \
--sample=WT281 --slide=V11F09-071 \
--slidefile=${slide_dir}/V11F09-071.gpr \
--area=B1 --localcores=2
echo 'WT281 complete'


echo 'Start 4WT'
spaceranger count --id=WT731 \
--transcriptome=${Ref} \
--fastqs=${fas_dir}/WT731 \
--image=${image_dir}/WT731.tif \
--sample=WT731 --slide=V11F09-071 \
--slidefile=${slide_dir}/V11F09-071.gpr \
--area=A1 --localcores=2
echo 'WT731 complete'
echo 'all complete'
