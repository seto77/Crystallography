"""空間群の回帰テスト用の基準データ (symmetry_reference.json) を生成する。

基準:
  - spglib  : Hall 番号 1..530 の対称操作 (センタリング込み) と空間群型の情報 (記号・点群)
  - cctbx   : Wyckoff 位置 (文字・多重度・サイト対称の型・特殊位置の射影演算子)。cctbx の操作は spglib と一致することを生成時に確認する
  - ReciPro 独自の通し番号 531..540 (センタリング付き三斜晶 C-1, I-1, A-1, B-1, F-1, A1, B1, C1, F1, I1) は cctbx で Hall 記号から生成する

ReciPro の通し番号 1..530 は spglib の Hall 番号と一致する (照合で確認済み)。

使い方:  pip install spglib cctbx-base
          python generate_symmetry_reference.py      (このファイルと同じフォルダに JSON を書き出す)
注意: cctbx は spglib より先に import すること (逆順だと環境によって segfault する)。
"""
import json, os, datetime
from cctbx import sgtbx
import spglib

DEN = 24  # 並進は 1/24 単位の整数で保存する (全設定で十分)
EXTRA = {531: '-C 1', 532: '-I 1', 533: '-A 1', 534: '-B 1', 535: '-F 1',
         536: 'A 1', 537: 'B 1', 538: 'C 1', 539: 'F 1', 540: 'I 1'}

def op_key(rot, tr):
    r = [int(round(v)) for v in rot]
    t = []
    for v in tr:
        n = round(v * DEN)
        if abs(v * DEN - n) > 1e-9:
            raise ValueError(f'translation {tr} is not a multiple of 1/{DEN}')
        t.append(int(n) % DEN)
    return r + t

def cctbx_ops(sg):
    out = []
    for op in sg.all_ops():
        out.append(op_key(op.r().as_double(), op.t().as_double()))
    return out

def wyckoff(sgi):
    tab = sgi.wyckoff_table()
    res = []
    for k in range(tab.size()):
        p = tab.position(k)
        so = p.special_op()
        res.append({'letter': p.letter(), 'mult': p.multiplicity(), 'site': str(p.point_group_type()),
                    'special': list(so.r().as_double()) + list(so.t().as_double())})
    return res

groups = []
for h in list(range(1, 531)) + sorted(EXTRA):
    if h <= 530:
        t = spglib.get_spacegroup_type(h)
        hall = t.hall_symbol
        s = spglib.get_symmetry_from_database(h)
        ops = [op_key(R.flatten(), tr) for R, tr in zip(s['rotations'], s['translations'])]
    else:
        hall = EXTRA[h]
    sgi = sgtbx.space_group_info(symbol='Hall: ' + hall)
    cops = cctbx_ops(sgi.group())
    if h > 530:
        ops = cops
    elif sorted(map(tuple, cops)) != sorted(map(tuple, ops)):
        raise RuntimeError(f'spglib and cctbx disagree for hall {h}')
    g = {'series': h, 'hall_symbol': hall, 'ops': sorted(ops),
         'cctbx_reference_setting': sgi.type().cb_op().is_identity_op(),
         'wyckoff': wyckoff(sgi)}
    if h <= 530:
        g.update({'number': t.number, 'hm_short': t.international_short, 'hm_full': t.international_full,
                  'schoenflies': t.schoenflies, 'choice': t.choice,
                  'pg_hm': t.pointgroup_international, 'pg_schoenflies': t.pointgroup_schoenflies})
    else:
        g.update({'number': 2 if hall.startswith('-') else 1})
    groups.append(g)

meta = {'generated': datetime.date.today().isoformat(), 'spglib': spglib.__version__,
        'cctbx': __import__('cctbx').__version__ if hasattr(__import__('cctbx'), '__version__') else 'unknown',
        'translation_denominator': DEN,
        'op_format': '[r11,r12,r13,r21,r22,r23,r31,r32,r33, t1*DEN, t2*DEN, t3*DEN] (x\' = R x + t, t mod 1)',
        'special_format': '[r11..r33, t1, t2, t3] projection operator onto the Wyckoff position (cctbx special_op)'}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'symmetry_reference.json')
# リポジトリの .editorconfig (json: utf-8-bom, crlf) に合わせる
with open(out, 'w', encoding='utf-8-sig', newline='') as f:
    json.dump({'meta': meta, 'groups': groups}, f, separators=(',', ':'))
    f.write('\r\n')
print(f'wrote {len(groups)} groups to {out}')
