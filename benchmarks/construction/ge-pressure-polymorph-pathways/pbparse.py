"""Minimal protobuf decoder for a TensorFlow GraphDef (no TensorFlow): {node name: (op, inputs, constant, attrs)}."""
import numpy as np, struct
def varint(b,i):
  r=0;s=0
  while True:
    c=b[i];i+=1;r|=(c&0x7f)<<s;s+=7
    if c<0x80: return r,i
def fields(b):
  i=0;out=[]
  while i<len(b):
    k,i=varint(b,i);f,t=k>>3,k&7
    if t==0: v,i=varint(b,i)
    elif t==1: v=b[i:i+8];i+=8
    elif t==5: v=b[i:i+4];i+=4
    elif t==2:
      l,i=varint(b,i);v=b[i:i+l];i+=l
    else: raise ValueError(t)
    out.append((f,t,v))
  return out
DT={1:np.float32,2:np.float64,3:np.int32,9:np.int64,7:None,10:np.bool_}
def tensor(b):
  dtype=None;shape=[];content=None;vals={}
  for f,t,v in fields(b):
    if f==1: dtype=v
    elif f==2:
      for f2,t2,v2 in fields(v):
        if f2==2:
          d=[x for x in fields(v2) if x[0]==1]; shape.append(d[0][2] if d else 0)
    elif f==4: content=v
    else: vals.setdefault(f,[]).append((t,v))
  if dtype==7:  # string
    return [v for t,v in vals.get(8,[])]
  npt=DT.get(dtype)
  if content is not None:
    a=np.frombuffer(content,dtype=npt)
  else:
    fld={1:5,2:6,3:7,9:10}.get(dtype)
    raw=vals.get(fld,[])
    xs=[]
    for t,v in raw:
      if t==2:
        if dtype in (1,): xs+=list(np.frombuffer(v,np.float32))
        elif dtype==2: xs+=list(np.frombuffer(v,np.float64))
        else:
          j=0
          while j<len(v):
            x,j=varint(v,j); xs.append(x-2**64 if x>=2**63 else x)
      elif t==5: xs.append(np.frombuffer(v,np.float32)[0])
      elif t==1: xs.append(np.frombuffer(v,np.float64)[0])
      elif t==0: xs.append(v-2**64 if v>=2**63 else v)
      else: xs.append(v)
    a=np.array(xs,dtype=npt)
  n=int(np.prod(shape)) if shape else 1
  if a.size==1 and n>1: a=np.full(n,a[0],dtype=npt)
  return a.reshape(shape) if a.size==n else a
def graph(path):
  return graph_from_bytes(open(path,'rb').read())
def graph_from_bytes(b):
  nodes={}
  for f,t,v in fields(b):
    if f!=1: continue
    name=op=None;inputs=[];attrs={}
    for f2,t2,v2 in fields(v):
      if f2==1: name=v2.decode()
      elif f2==2: op=v2.decode()
      elif f2==3: inputs.append(v2.decode())
      elif f2==5:
        kv=fields(v2); key=kv[0][2].decode(); val=kv[1][2] if len(kv)>1 else b''
        attrs[key]=val
    const=None
    if op=='Const' and 'value' in attrs:
      for f3,t3,v3 in fields(attrs['value']):
        if f3==8: const=tensor(v3)
    nodes[name]=(op,inputs,const,attrs)
  return nodes
