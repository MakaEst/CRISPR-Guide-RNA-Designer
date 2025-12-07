"""
CRISPR gRNA Designer — простая реализация для SpCas9 (NGG PAM).
Поддерживает поиск PAM на обеих цепях, генерацию guide, базовый on-target scoring
и примитивный off-target поиск по локальному референсу (FASTA) методом скользящего окна.
"""


from typing import List, Dict, Tuple
import re
from Bio import SeqIO
import pandas as pd
import numpy as np


# Небольшкие утилиты
BASES = set('ACGT')


def revcomp(seq: str) -> str:
comp = str.maketrans('ACGTacgt', 'TGCAtgca')
return seq.translate(comp)[::-1]




def _pam_regex(pam_pattern: str) -> str:
# Только поддержка N, A,C,G,T и конкретных букв: пример NGG -> .GG
pat = pam_pattern.upper().replace('N', '.')
return pat




def find_guides(seq: str, pam_pattern: str='NGG', guide_len: int=20) -> pd.DataFrame:
"""Ищет PAM на прямой и обратной цепях, возвращает кандидаты gRNA.
Возвращаемые поля: guide_seq, pam, start, end, strand
"""
s = seq.upper()
pam_re = _pam_regex(pam_pattern)
results = []
# Прямая цепь: guide находится непосредственно перед PAM (5'->3')
for m in re.finditer(pam_re, s):
pam_start = m.start()
guide_start = pam_start - guide_len
if guide_start >= 0:
guide = s[guide_start:guide_start+guide_len]
if set(guide).issubset(BASES):
results.append({
'guide_seq': guide,
'pam': s[pam_start:pam_start+len(m.group(0))],
'start': guide_start,
'end': pam_start+len(m.group(0)),
'strand': '+'
})
# Обратная цепь: ищем PAM на обратной цепи
rs = revcomp(s)
for m in re.finditer(pam_re, rs):
pam_start = m.start()
guide_start = pam_start - guide_len
if guide_start >= 0:
guide = rs[guide_start:guide_start+guide_len]
if set(guide).issubset(BASES):
# Переводим координаты обратно в оригинальную последовательность
# в оригинальной (+) координате: start = len(s) - (pam_start+len(pam))
pam_len = len(m.group(0))
orig_pam_end = len(s) - pam_start
orig_guide_start = len(s) - (pam_start+pam_len) - guide_len
results.append({
'guide_seq': guide,
'pam': rs[pam_start:pam_start+pam_len],
'start': orig_guide_start,
'end': orig_pam_end,
'strand': '-'
})
df = pd.DataFrame(results)
print(scored[['guide_seq','pam','start','end','strand','gc_pct','on_target_score','off_target_penalty','final_score']].to_string(index=False))