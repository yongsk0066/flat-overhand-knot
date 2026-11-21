# 논문용: 느슨한 유효 4_1 접힘 상태 (조임 전) — 영문, 무제목 — for the paper: a loose valid 4_1 folded state (before tightening) — English, untitled
import sys, json
sys.path.insert(0,'verify')
from n4_render import render
h=json.load(open("verify/n4_valid41.json"))
th,o=h[1]
json.dump({"theta":th,"order":o},open("/tmp/loose_state.json","w"))
