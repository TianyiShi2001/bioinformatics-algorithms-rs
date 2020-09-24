use std::rc::Rc;

fn main() {
    println!("{:?}", "hello");
}

enum Operation {
    Insert,
    Delete,
    Match,
    Mismatch,
}

enum TbCell {
    Cons(Vec<(Operation, Rc<TbCell>)>),
    Nil,
}

impl TbCell {
    pub fn traceback(&self) {
        let l = TbCell::Nil;
    }

    pub fn traceback_recursive(curr: &TbCell, res: TbCell) -> TbCell {
        match curr {
            TbCell::Cons(v) => {
                let w = Vec::new();
                for (operation, tbcell) in v {}
            }
            TbCell::Nil => return res,
        }
    }
}
